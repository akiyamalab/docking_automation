"""
N×M ドッキングPoC スクリプト

N個のタンパク質 × M個の化合物のドッキング計算を実行し、
メトリクスをJSON Lines形式で記録する。HDF5リポジトリによる冪等性を保証。

Usage:
    cd /workspaces/20260422_mouse_docking/docking_automation
    python3 examples/nxm_poc.py
    python3 examples/nxm_poc.py --dry-run
    python3 examples/nxm_poc.py --max-compounds 5
"""

import argparse
import gzip
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Optional

import numpy as np
import psutil

os.environ["BABEL_QUIET"] = "1"

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent

INPUT_DIR = SCRIPT_DIR / "input"
OUTPUT_DIR = SCRIPT_DIR / "output"
ALPHAFOLD_DIR = INPUT_DIR / "alphafold"
AFDB_MOUSE_DIR = INPUT_DIR / "afdb_mouse"
ALDR_DIR = INPUT_DIR / "ALDR"

ACTIVES_SUBSET_SDF = ALDR_DIR / "actives_subset.sdf"
ACTIVES_FINAL_GZ = ALDR_DIR / "actives_final.sdf.gz"
ACTIVES_EXTENDED20_SDF = ALDR_DIR / "actives_extended20.sdf"
PROTEIN_LIST_JSON = AFDB_MOUSE_DIR / "protein_list.json"

METRICS_JSONL = OUTPUT_DIR / "nxm_poc_metrics.jsonl"
SUMMARY_TXT = OUTPUT_DIR / "nxm_poc_summary.txt"
HDF5_DIR = OUTPUT_DIR / "nxm_poc_hdf5"


def grid_box_centroid(pdb_path: Path, size: float = 30.0):
    """タンパク質の重心を中心とした固定サイズのグリッドボックスを生成する。

    PDBファイルのCα原子座標から重心を計算する。
    """
    from docking_automation.docking.grid_box import GridBox

    coords = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue

    if not coords:
        # CA原子が見つからない場合は全ATOM座標を使う
        with open(pdb_path) as f:
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append([x, y, z])
                except ValueError:
                    continue

    if not coords:
        raise ValueError(f"座標が見つかりません: {pdb_path}")

    center = np.mean(np.array(coords), axis=0)
    return GridBox(center=tuple(center), size=(size, size, size))


def prepare_extended20() -> Path:
    """actives_final.sdf.gz から先頭20件を抽出して actives_extended20.sdf を生成する。"""
    from rdkit import Chem
    from rdkit.Chem import SDWriter

    if ACTIVES_EXTENDED20_SDF.exists():
        return ACTIVES_EXTENDED20_SDF

    print("actives_extended20.sdf を生成中...")
    writer = SDWriter(str(ACTIVES_EXTENDED20_SDF))
    count = 0
    with gzip.open(ACTIVES_FINAL_GZ, "rb") as f:
        suppl = Chem.ForwardSDMolSupplier(f)
        for mol in suppl:
            if mol is None:
                continue
            writer.write(mol)
            count += 1
            if count >= 20:
                break
    writer.close()
    print(f"actives_extended20.sdf 生成完了: {count}件")
    return ACTIVES_EXTENDED20_SDF


def _count_residues(pdb_path: Path) -> int:
    """BioPythonでPDBファイルの残基数を返す。"""
    from Bio import PDB
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", str(pdb_path))
    return sum(1 for _ in structure.get_residues())


def load_proteins(protein_list_json: Optional[Path] = None) -> List:
    """protein_list.json からタンパク質リストを読み込む。

    ファイルが存在しない場合は alphafold/ ディレクトリの2件にフォールバックする。
    残基数 > 500 のタンパク質は警告ログを出してスキップする。
    """
    from docking_automation.molecule.protein import Protein

    proteins = []
    json_path = protein_list_json or PROTEIN_LIST_JSON

    if json_path.exists():
        with open(json_path) as f:
            paths = json.load(f)
        for p in paths:
            path = Path(p)
            if path.exists():
                residue_count = _count_residues(path)
                if residue_count > 500:
                    print(f"[WARN] 残基数 {residue_count} > 500、スキップ: {path.name}")
                    continue
                protein_id = path.stem
                # AF-XXXXX-F1-model_v4 → AF-XXXXX
                if "-F1-" in protein_id:
                    protein_id = protein_id.split("-F1-")[0]
                proteins.append(Protein(path, id=protein_id))
            else:
                print(f"[WARN] タンパク質ファイルが見つかりません: {path}")
    else:
        print(f"[WARN] {json_path} が見つかりません。alphafold/ ディレクトリにフォールバックします。")
        for pdb in sorted(ALPHAFOLD_DIR.glob("*.pdb")):
            residue_count = _count_residues(pdb)
            if residue_count > 500:
                print(f"[WARN] 残基数 {residue_count} > 500、スキップ: {pdb.name}")
                continue
            protein_id = pdb.stem
            if "-F1-" in protein_id:
                protein_id = protein_id.split("-F1-")[0]
            proteins.append(Protein(pdb, id=protein_id))

    print(f"タンパク質: {len(proteins)}件 読み込み完了")
    return proteins


def load_compounds(sdf_path: Path, max_n: Optional[int] = 10):
    """SDF ファイルから化合物セットを読み込む。

    max_n が指定された場合は先頭 max_n 件に制限する。
    """
    from docking_automation.molecule.compound_set import CompoundSet

    compound_set = CompoundSet(sdf_path)
    total = compound_set.get_compound_count()
    print(f"化合物: {total}件 (ファイル: {sdf_path.name})")

    if max_n is not None and max_n < total:
        indices = set(range(max_n))
        compound_set = compound_set.with_indices(indices)
        print(f"先頭 {max_n} 件に制限")

    return compound_set


def run_nxm(proteins: List, compound_set, repo, metrics_fp) -> List[dict]:
    """N×M ドッキングを実行してメトリクスを記録する。"""
    from vina import Vina

    from docking_automation.converters.molecule_converter import MoleculeConverter
    from docking_automation.docking.autodockvina_docking import AutoDockVina, AutoDockVinaParameters
    from docking_automation.docking.docking_result import DockingResult

    docking_tool = AutoDockVina()
    converter = MoleculeConverter()
    process = psutil.Process()
    all_metrics = []

    try:
        from tqdm import tqdm
        protein_iter = tqdm(proteins, desc="Proteins")
    except ImportError:
        protein_iter = proteins

    for protein in protein_iter:
        protein_id = protein.id
        print(f"\n[Protein] {protein_id}")

        # --- タンパク質前処理 ---
        t_pre_start = time.time()
        try:
            prep_protein = docking_tool._preprocess_protein(protein)
        except Exception as e:
            print(f"  [ERROR] タンパク質前処理失敗: {e}")
            continue
        preprocess_sec = time.time() - t_pre_start

        # --- グリッドボックス計算 ---
        t_grid_start = time.time()
        try:
            grid_box = grid_box_centroid(protein.path)
        except Exception as e:
            print(f"  [ERROR] グリッドボックス計算失敗: {e}")
            continue
        grid_sec = time.time() - t_grid_start

        center = grid_box.center
        size = grid_box.size

        # --- 化合物前処理 ---
        try:
            prep_compounds = docking_tool._preprocess_compound_set(compound_set)
        except Exception as e:
            print(f"  [ERROR] 化合物前処理失敗: {e}")
            continue

        try:
            from tqdm import tqdm
            compound_iter = tqdm(
                enumerate(prep_compounds.file_paths),
                total=len(prep_compounds.file_paths),
                desc="  Compounds",
                leave=False,
            )
        except ImportError:
            compound_iter = enumerate(prep_compounds.file_paths)

        for idx, compound_path in compound_iter:
            try:
                compound_hash = prep_compounds.get_compound_hash(idx)
                mem_mb = process.memory_info().rss / 1024 / 1024
                t_dock_start = time.time()

                if repo._exists(protein.content_hash, compound_hash):
                    result = repo.load_by_hashes(protein.content_hash, compound_hash)
                    dock_sec = time.time() - t_dock_start
                    hdf5_write_sec = 0.0
                    reused = True
                    score = result.docking_score if result is not None else None
                else:
                    # Vina でドッキング実行
                    temp_dir = Path(tempfile.mkdtemp())
                    output_pdbqt = temp_dir / f"output_{idx}.pdbqt"
                    output_sdf = temp_dir / f"output_{idx}.sdf"

                    v = Vina(cpu=1, seed=1, verbosity=0)
                    v.set_receptor(str(prep_protein.file_path))
                    v.set_ligand_from_file(str(compound_path))
                    v.compute_vina_maps(
                        center=[float(center[0]), float(center[1]), float(center[2])],
                        box_size=[float(size[0]), float(size[1]), float(size[2])],
                    )
                    v.dock(exhaustiveness=4, n_poses=3, min_rmsd=1.0)
                    v.write_poses(str(output_pdbqt), n_poses=3, overwrite=True)
                    converter.pdbqt_to_sdf(output_pdbqt, output_sdf)
                    scores = v.energies()
                    dock_sec = time.time() - t_dock_start

                    result = DockingResult(
                        result_path=output_sdf,
                        protein_id=protein_id,
                        compound_set_id=compound_path.stem.split("_")[0],
                        compound_index=idx,
                        docking_score=float(scores[0, 0]),
                        protein_content_hash=protein.content_hash,
                        compound_content_hash=compound_hash,
                        compoundset_content_hash=prep_compounds.content_hash,
                        metadata={"tool": "AutoDock Vina", "exhaustiveness": 4, "num_modes": 3},
                    )

                    # HDF5 保存
                    t_write_start = time.time()
                    repo.save(result)
                    hdf5_write_sec = time.time() - t_write_start
                    reused = False
                    score = float(scores[0, 0])

                metric = {
                    "protein_id": protein_id,
                    "compound_idx": idx,
                    "preprocess_sec": round(preprocess_sec, 3),
                    "grid_sec": round(grid_sec, 4),
                    "dock_sec": round(dock_sec, 3),
                    "hdf5_write_sec": round(hdf5_write_sec, 4),
                    "memory_mb": round(mem_mb, 1),
                    "score": score,
                    "reused": reused,
                    "timestamp": int(time.time()),
                }
                metrics_fp.write(json.dumps(metric, ensure_ascii=False) + "\n")
                metrics_fp.flush()
                all_metrics.append(metric)

                status = "reused" if reused else f"score={score:.2f}"
                print(f"    compound[{idx}]: {status} ({dock_sec:.1f}s)")

            except Exception as e:
                print(f"  [ERROR] compound[{idx}] 処理失敗: {e}")
                continue

    return all_metrics


def write_summary(all_metrics: List[dict], elapsed_total: float) -> None:
    """サマリファイルを書き出す。"""
    if not all_metrics:
        with open(SUMMARY_TXT, "w") as f:
            f.write("メトリクスなし（エラーで全ペア失敗）\n")
        return

    n_total = len(all_metrics)
    n_reused = sum(1 for m in all_metrics if m["reused"])
    n_new = n_total - n_reused
    scores = [m["score"] for m in all_metrics if m["score"] is not None]
    avg_score = np.mean(scores) if scores else None
    avg_dock = np.mean([m["dock_sec"] for m in all_metrics])
    avg_preprocess = np.mean([m["preprocess_sec"] for m in all_metrics])
    bottleneck = "前処理" if avg_preprocess > avg_dock else "ドッキング"

    lines = [
        "=== N×M ドッキングPoC サマリ ===",
        f"総ペア数:        {n_total}",
        f"新規ドッキング:  {n_new}",
        f"再利用:          {n_reused}",
        f"総経過時間:      {elapsed_total:.1f}s",
        f"平均ドッキング時間: {avg_dock:.2f}s/ペア",
        f"平均前処理時間:  {avg_preprocess:.2f}s/タンパク質",
        f"平均スコア:      {avg_score:.3f}" if avg_score is not None else "平均スコア: N/A",
        f"ボトルネック:    {bottleneck}",
    ]
    text = "\n".join(lines) + "\n"
    print("\n" + text)
    with open(SUMMARY_TXT, "w") as f:
        f.write(text)


def parse_args():
    parser = argparse.ArgumentParser(description="N×M ドッキングPoC")
    parser.add_argument("--dry-run", action="store_true", help="設定確認のみ（ドッキング実行なし）")
    parser.add_argument("--max-compounds", type=int, default=10, help="処理する化合物の最大数（デフォルト: 10）")
    parser.add_argument("--extended", action="store_true", help="actives_extended20.sdf を使用（20件）")
    parser.add_argument(
        "--protein-list", type=Path, default=None, help="protein_list.json のパス（デフォルト: afdb_mouse/protein_list.json）"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== N×M ドッキングPoC 開始 ===")
    print(f"作業ディレクトリ: {PROJECT_DIR}")

    # --- actives_extended20.sdf 準備 ---
    if ACTIVES_FINAL_GZ.exists():
        prepare_extended20()

    # --- 入力ファイル選択 ---
    if args.extended:
        compound_sdf = ACTIVES_EXTENDED20_SDF
        max_n = args.max_compounds if args.max_compounds != 10 else 20
    else:
        compound_sdf = ACTIVES_SUBSET_SDF
        max_n = args.max_compounds

    if not compound_sdf.exists():
        print(f"[ERROR] 化合物ファイルが見つかりません: {compound_sdf}")
        sys.exit(1)

    # --- タンパク質と化合物の読み込み ---
    proteins = load_proteins(args.protein_list)
    if not proteins:
        print("[ERROR] タンパク質が1件も読み込めませんでした。")
        sys.exit(1)

    compound_set = load_compounds(compound_sdf, max_n=max_n)

    print(f"\n実行設定:")
    print(f"  タンパク質数:   {len(proteins)}")
    print(f"  化合物数:       {max_n}")
    print(f"  総ペア数:       {len(proteins) * max_n}")
    print(f"  exhaustiveness: 4")
    print(f"  num_modes:      3")
    print(f"  HDF5 出力先:    {HDF5_DIR}")

    if args.dry_run:
        print("\n[dry-run] 設定確認完了。ドッキングは実行しません。")
        return

    # --- 出力ディレクトリ準備 ---
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    HDF5_DIR.mkdir(parents=True, exist_ok=True)

    # --- HDF5 リポジトリ初期化 ---
    from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
        DockingResultRepositoryFactory,
        RepositoryType,
    )

    repo = DockingResultRepositoryFactory.create(
        RepositoryType.HDF5,
        base_directory=HDF5_DIR,
        config={"mode": "append"},
    )

    # --- N×M ドッキング実行 ---
    t_start = time.time()
    with open(METRICS_JSONL, "w") as metrics_fp:
        all_metrics = run_nxm(proteins, compound_set, repo, metrics_fp)
    elapsed = time.time() - t_start

    # --- サマリ出力 ---
    write_summary(all_metrics, elapsed)
    print(f"\nメトリクス: {METRICS_JSONL}")
    print(f"サマリ:     {SUMMARY_TXT}")
    print("=== 完了 ===")


if __name__ == "__main__":
    main()
