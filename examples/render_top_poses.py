"""上位ドッキングポーズの PyMOL バッチ描画。

入力ソース:
  1. HDF5 リポジトリ — pose_blob (gzip SDF) を抽出 → PyMOL
  2. Uni-Dock 生出力ディレクトリ (*_out.pdbqt) — PDBQT を直接 PyMOL に読ませる

使用例 (headless レンダリング):
    /opt/miniforge/bin/python examples/render_top_poses.py \
        --pdbqt-root /tmp/mps_bench/grid/N4_L160/run_0 \
        --receptor /tmp/mps_bench/receptors/AF-A0A0A2JW93-F1-model_v4.pdbqt \
        --top-k 10

PyMOL は headless 利用。X11 不要。`ray` によるオフスクリーン描画で PNG 出力。
"""
from __future__ import annotations

import argparse
import gzip
import re
import tempfile
from pathlib import Path
from typing import List, Optional, Tuple

VINA_RESULT_RE = re.compile(r"REMARK VINA RESULT:\s+(-?\d+\.\d+)")


def parse_pdbqt_score(pdbqt: Path) -> Optional[float]:
    try:
        with pdbqt.open() as f:
            for line in f:
                m = VINA_RESULT_RE.match(line)
                if m:
                    return float(m.group(1))
    except OSError:
        return None
    return None


def collect_top_from_pdbqt(root: Path, k: int,
                           lo: float = -1e9, hi: float = 1e9) -> List[Tuple[float, Path]]:
    scored: List[Tuple[float, Path]] = []
    for p in root.rglob("*_out.pdbqt"):
        s = parse_pdbqt_score(p)
        if s is not None and lo <= s <= hi:
            scored.append((s, p))
    scored.sort(key=lambda x: x[0])
    return scored[:k]


def collect_top_from_hdf5(hdf5_dir: Path, k: int) -> List[Tuple[float, Path, str, str]]:
    """Return list of (score, sdf_tempfile, protein_id, compound_hash)."""
    import h5py

    buf: List[Tuple[float, bytes, str, str]] = []
    for h5 in sorted(hdf5_dir.glob("*.h5")):
        with h5py.File(h5, "r") as f:
            if "results" not in f:
                continue
            for ph in f["results"]:
                for ch in f[f"results/{ph}"]:
                    g = f[f"results/{ph}/{ch}"]
                    score = float(g["docking_score"][()])
                    if len(buf) >= k and score >= buf[-1][0]:
                        continue
                    pose_bytes = gzip.decompress(bytes(g["pose_blob"][()]))
                    pid = g.attrs.get("protein_id", "")
                    buf.append((score, pose_bytes, pid, ch))
                    buf.sort(key=lambda x: x[0])
                    buf = buf[:k]

    out: List[Tuple[float, Path, str, str]] = []
    for score, pose_bytes, pid, ch in buf:
        tmp = Path(tempfile.mkstemp(suffix=".sdf")[1])
        tmp.write_bytes(pose_bytes)
        out.append((score, tmp, pid, ch))
    return out


def render_pose(receptor: Path, ligand: Path, score: float, out_png: Path, label: str) -> None:
    from pymol import cmd

    cmd.reinitialize()
    cmd.load(str(receptor), "rec")
    cmd.load(str(ligand), "lig")

    cmd.hide("everything", "rec")
    cmd.show("cartoon", "rec")
    cmd.color("gray80", "rec")
    cmd.show("sticks", "lig")
    cmd.color("cyan", "lig and elem C")
    cmd.show("surface", "rec within 5 of lig")
    cmd.set("transparency", 0.6)

    cmd.zoom("lig", 6)
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 1)
    cmd.set("ray_shadow", 0)

    cmd.label("lig and name C1", f'"{label} ({score:.2f})"')
    cmd.ray(900, 700)
    cmd.png(str(out_png), dpi=150)


def main():
    ap = argparse.ArgumentParser()
    src = ap.add_mutually_exclusive_group(required=True)
    src.add_argument("--hdf5-dir", type=Path)
    src.add_argument("--pdbqt-root", type=Path)
    ap.add_argument("--receptor", type=Path,
                    help="(PDBQT root 使用時) 単一の receptor PDBQT/PDB")
    ap.add_argument("--receptor-dir", type=Path,
                    help="(HDF5 使用時) protein_id から receptor PDB を引くディレクトリ")
    ap.add_argument("--top-k", type=int, default=10)
    ap.add_argument("--min-score", type=float, default=-20.0)
    ap.add_argument("--max-score", type=float, default=5.0)
    ap.add_argument("--out-dir", type=Path,
                    default=Path(__file__).parent / "output" / "analysis" / "poses")
    args = ap.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.pdbqt_root:
        if args.receptor is None:
            ap.error("--pdbqt-root 使用時は --receptor も必須")
        top = collect_top_from_pdbqt(args.pdbqt_root, args.top_k,
                                     args.min_score, args.max_score)
        print(f"PDBQT ツリーから top-{args.top_k}:")
        for i, (score, p) in enumerate(top):
            label = p.stem
            out_png = args.out_dir / f"top_{i:02d}_{label}.png"
            print(f"  [{i:2d}] {label}  score={score:.2f}  -> {out_png.name}")
            render_pose(args.receptor, p, score, out_png, label)
    else:
        if args.receptor_dir is None:
            ap.error("--hdf5-dir 使用時は --receptor-dir も必須")
        top = collect_top_from_hdf5(args.hdf5_dir, args.top_k)
        print(f"HDF5 から top-{args.top_k}:")
        for i, (score, sdf, pid, ch) in enumerate(top):
            rec_candidates = list(args.receptor_dir.glob(f"*{pid}*.pdb"))
            if not rec_candidates:
                print(f"  [skip] receptor for {pid} not found")
                continue
            out_png = args.out_dir / f"top_{i:02d}_{pid}_{ch[:8]}.png"
            label = f"{pid}/{ch[:8]}"
            print(f"  [{i:2d}] {label}  score={score:.2f}")
            render_pose(rec_candidates[0], sdf, score, out_png, label)

    print(f"\n描画完了: {args.out_dir}")


if __name__ == "__main__":
    main()
