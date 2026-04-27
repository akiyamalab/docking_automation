#!/usr/bin/env python3
"""HDF5 ドッキング結果から指定範囲を「receptor PDB + multi-record pose SDF」に書き出す。

出力レイアウト:
    <out>/
      <protein_id>/
        receptor.pdb     # --receptors-dir から該当ファイルをコピー
        poses.sdf        # 該当 protein × 選択 compound の pose を全件束ねた multi-record SDF
                         # 各 record に protein_id / compound_hash / compound_set_id /
                         # compound_index / docking_score を property として付与

範囲指定 (AND で適用):
    --proteins / --proteins-file       protein_id (HDF5 attrs.protein_id) を絞り込む
    --compounds-file                   "compound_set_id<TAB>compound_index" or "compound_hash" の行
    --top-k                            per-protein で score 昇順 top-K を取得
    --score-max                        score <= X の pose のみ
    --limit-pairs                      安全上限 (default 1000)。超過時は早期 abort

対応スキーマ: v2 (per-pair: /results/<phash>/<chash>/{pose_blob,docking_score}),
              v3 (protein-bundle: /results/<phash>/{compound_hashes[],pose_blobs[],docking_scores[]})

Usage:
    python scripts/export_poses.py \\
        --hdf5 results/run.h5 \\
        --receptors-dir jobs/inputs/receptors_pdb \\
        --proteins AF-XXXX-F1,AF-YYYY-F1 \\
        --top-k 5 \\
        --out exported/
"""
from __future__ import annotations

import argparse
import gzip
import logging
import shutil
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

import h5py

logger = logging.getLogger("export_poses")

_FAILED_SCORE_SENTINEL = -999.0


def _read_lines(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip() and not ln.startswith("#")]


def _parse_compound_filter(path: Optional[Path]) -> Tuple[Set[str], Set[Tuple[str, int]]]:
    """compounds-file を読み込む。

    各行は以下のいずれか:
      - "compound_hash"           (single token, 任意長 hash)
      - "compound_set_id<TAB>compound_index"
    """
    hashes: Set[str] = set()
    set_idx: Set[Tuple[str, int]] = set()
    if path is None:
        return hashes, set_idx
    for line in _read_lines(path):
        parts = line.split("\t") if "\t" in line else line.split()
        if len(parts) == 1:
            hashes.add(parts[0])
        elif len(parts) >= 2:
            try:
                set_idx.add((parts[0], int(parts[1])))
            except ValueError:
                logger.warning(f"compound 行を解析できません: {line!r}")
    return hashes, set_idx


def _decompress_pose(blob: bytes) -> Optional[str]:
    if not blob:
        return None
    try:
        return gzip.decompress(blob).decode("utf-8")
    except OSError:
        return None


def _iter_v2(f: h5py.File) -> Iterable[Dict[str, Any]]:
    if "results" not in f:
        return
    for ph in f["results"]:
        pgrp = f[f"results/{ph}"]
        for ch in pgrp:
            g = pgrp[ch]
            if "docking_score" not in g:
                continue
            score = float(g["docking_score"][()])
            blob = bytes(g["pose_blob"][()]) if "pose_blob" in g else b""
            yield {
                "protein_hash": ph,
                "compound_hash": ch,
                "score": score,
                "pose_blob": blob,
                "protein_id": g.attrs.get("protein_id", ""),
                "compound_set_id": g.attrs.get("compound_set_id", ""),
                "compound_index": int(g.attrs.get("compound_index", -1)),
            }


def _iter_v3(f: h5py.File) -> Iterable[Dict[str, Any]]:
    if "results" not in f:
        return
    for ph in f["results"]:
        grp = f[f"results/{ph}"]
        if "compound_hashes" not in grp:
            continue
        n = grp["compound_hashes"].shape[0]
        hashes = grp["compound_hashes"][:]
        scores = grp["docking_scores"][:]
        for i in range(n):
            score = float(scores[i])
            if score == _FAILED_SCORE_SENTINEL:
                continue
            ch = hashes[i].decode("utf-8") if isinstance(hashes[i], bytes) else str(hashes[i])
            yield {
                "protein_hash": ph,
                "compound_hash": ch,
                "score": score,
                "pose_blob": bytes(grp["pose_blobs"][i]),
                "protein_id": "",
                "compound_set_id": "",
                "compound_index": -1,
            }


def _detect_schema(f: h5py.File) -> str:
    if "results" not in f:
        return "empty"
    for ph in f["results"]:
        grp = f[f"results/{ph}"]
        if "compound_hashes" in grp:
            return "v3"
        for ch in grp:
            return "v2"
    return "empty"


def _apply_filters(
    entries: Iterable[Dict[str, Any]],
    proteins: Optional[Set[str]],
    compound_hashes: Set[str],
    set_idx: Set[Tuple[str, int]],
    score_max: Optional[float],
) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    for e in entries:
        if proteins is not None and e["protein_id"] and e["protein_id"] not in proteins:
            continue
        if proteins is not None and not e["protein_id"]:
            # v3 では attrs.protein_id がないので protein_hash でも一致を許す
            if e["protein_hash"] not in proteins:
                continue
        if compound_hashes and e["compound_hash"] not in compound_hashes:
            continue
        if set_idx and (e["compound_set_id"], e["compound_index"]) not in set_idx:
            continue
        if score_max is not None and e["score"] > score_max:
            continue
        out.append(e)
    return out


def _per_protein_top_k(entries: List[Dict[str, Any]], k: Optional[int]) -> List[Dict[str, Any]]:
    if k is None:
        return entries
    by_protein: Dict[str, List[Dict[str, Any]]] = {}
    for e in entries:
        key = e["protein_id"] or e["protein_hash"]
        by_protein.setdefault(key, []).append(e)
    out: List[Dict[str, Any]] = []
    for v in by_protein.values():
        v.sort(key=lambda x: x["score"])
        out.extend(v[:k])
    return out


def _copy_receptor(protein_key: str, receptors_dir: Path, out_dir: Path) -> Optional[Path]:
    """receptors_dir から protein_key にマッチする PDB を探してコピー。

    マッチ規則 (順に試行):
      1. <protein_key>.pdb
      2. *<protein_key>*.pdb (1 件のみヒットすれば採用)
    """
    cand = receptors_dir / f"{protein_key}.pdb"
    if cand.exists():
        target = out_dir / "receptor.pdb"
        shutil.copy2(cand, target)
        return target
    matches = list(receptors_dir.glob(f"*{protein_key}*.pdb"))
    if len(matches) == 1:
        target = out_dir / "receptor.pdb"
        shutil.copy2(matches[0], target)
        return target
    if len(matches) > 1:
        logger.warning(f"{protein_key}: receptor PDB が複数 hit ({len(matches)} 件), スキップ")
    else:
        logger.warning(f"{protein_key}: receptor PDB が {receptors_dir} に見つかりません")
    return None


def _write_protein_sdf(entries: List[Dict[str, Any]], out_path: Path) -> int:
    """RDKit を使い、protein 単位の multi-record SDF を書き出す。

    Returns: 実際に書き込まれた pose 数 (RDKit が parse 失敗した record はスキップ)。
    """
    from rdkit import Chem

    written = 0
    with Chem.SDWriter(str(out_path)) as w:
        for e in sorted(entries, key=lambda x: x["score"]):
            sdf_text = _decompress_pose(e["pose_blob"])
            if not sdf_text:
                continue
            suppl = Chem.SDMolSupplier()
            suppl.SetData(sdf_text, removeHs=False)
            for mol in suppl:
                if mol is None:
                    continue
                mol.SetProp("docking_score", f"{e['score']:.3f}")
                if e["protein_id"]:
                    mol.SetProp("protein_id", e["protein_id"])
                mol.SetProp("protein_hash", e["protein_hash"])
                mol.SetProp("compound_hash", e["compound_hash"])
                if e["compound_set_id"]:
                    mol.SetProp("compound_set_id", e["compound_set_id"])
                if e["compound_index"] >= 0:
                    mol.SetProp("compound_index", str(e["compound_index"]))
                w.write(mol)
                written += 1
    return written


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hdf5", type=Path, required=True, help="入力 HDF5")
    ap.add_argument("--receptors-dir", type=Path, required=True, help="受容体 PDB の置き場")
    ap.add_argument("--out", type=Path, required=True, help="出力ディレクトリ")
    ap.add_argument("--proteins", type=str, default=None, help="protein_id をカンマ区切りで指定")
    ap.add_argument("--proteins-file", type=Path, default=None, help="protein_id を 1 行 1 件で記載したファイル")
    ap.add_argument("--compounds-file", type=Path, default=None,
                    help="compound 絞り込み: 'compound_hash' or 'compound_set_id<TAB>compound_index'")
    ap.add_argument("--top-k", type=int, default=None, help="per-protein で score 昇順 top-K")
    ap.add_argument("--score-max", type=float, default=None, help="score <= X のみ")
    ap.add_argument("--limit-pairs", type=int, default=1000,
                    help="安全上限 (default 1000)。超過時は abort")
    ap.add_argument("--schema", choices=["v2", "v3", "auto"], default="auto")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args(argv)

    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO,
                        format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    if not args.hdf5.exists():
        logger.error(f"HDF5 が存在しません: {args.hdf5}")
        return 2
    if not args.receptors_dir.is_dir():
        logger.error(f"receptors-dir が存在しません: {args.receptors_dir}")
        return 2

    proteins: Optional[Set[str]] = None
    if args.proteins:
        proteins = set(p.strip() for p in args.proteins.split(",") if p.strip())
    if args.proteins_file:
        proteins = (proteins or set()) | set(_read_lines(args.proteins_file))

    compound_hashes, set_idx = _parse_compound_filter(args.compounds_file)

    with h5py.File(args.hdf5, "r") as f:
        schema = args.schema if args.schema != "auto" else _detect_schema(f)
        if schema == "empty":
            logger.error("HDF5 に results グループがありません")
            return 2
        logger.info(f"schema = {schema}")
        if schema == "v2":
            entries = list(_iter_v2(f))
        else:
            entries = list(_iter_v3(f))

    logger.info(f"全 entry: {len(entries)}")
    entries = _apply_filters(entries, proteins, compound_hashes, set_idx, args.score_max)
    logger.info(f"フィルタ後: {len(entries)}")
    entries = _per_protein_top_k(entries, args.top_k)
    logger.info(f"top-k 後: {len(entries)}")

    if len(entries) > args.limit_pairs:
        logger.error(
            f"対象 pose 数 {len(entries)} が --limit-pairs={args.limit_pairs} を超えています。"
            f" --top-k / --score-max / --proteins 等で絞り込むか、--limit-pairs を増やしてください。"
        )
        return 3
    if not entries:
        logger.warning("該当 pose がありません。条件を見直してください。")
        return 0

    by_protein: Dict[str, List[Dict[str, Any]]] = {}
    for e in entries:
        by_protein.setdefault(e["protein_id"] or e["protein_hash"], []).append(e)

    args.out.mkdir(parents=True, exist_ok=True)
    total_pose = 0
    skipped_no_pdb = 0
    for pkey, es in sorted(by_protein.items()):
        pdir = args.out / pkey
        pdir.mkdir(parents=True, exist_ok=True)
        rec = _copy_receptor(pkey, args.receptors_dir, pdir)
        if rec is None:
            skipped_no_pdb += 1
        n = _write_protein_sdf(es, pdir / "poses.sdf")
        total_pose += n
        logger.info(f"{pkey}: {n} pose 書き出し → {pdir}")

    logger.info(f"完了: protein={len(by_protein)}, pose={total_pose}, receptor 不明={skipped_no_pdb}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
