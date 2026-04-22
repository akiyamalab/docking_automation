"""
N×M ドッキング並列化PoC スクリプト (Dask LocalCluster版)

Dask LocalCluster によるタンパク質レベル並列化。
1タンパク質 × 全化合物 = 1 Daskタスク単位。
ワーカーはHDF5を開かない（結果収集パターン）。

Usage:
    cd /workspaces/20260422_mouse_docking/docking_automation
    python3 examples/nxm_poc_parallel.py --dry-run
    python3 examples/nxm_poc_parallel.py --workers 2 --exhaustiveness 1 --max-compounds 2
    python3 examples/nxm_poc_parallel.py --workers 4 --exhaustiveness 1
"""

import argparse
import json
import os
import sys
import tempfile
import time
from pathlib import Path
from typing import List, Optional

import numpy as np

os.environ["BABEL_QUIET"] = "1"

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent

sys.path.insert(0, str(SCRIPT_DIR))

INPUT_DIR = SCRIPT_DIR / "input"
OUTPUT_DIR = SCRIPT_DIR / "output"
AFDB_MOUSE_DIR = INPUT_DIR / "afdb_mouse"
ALDR_DIR = INPUT_DIR / "ALDR"

ACTIVES_SUBSET_SDF = ALDR_DIR / "actives_subset.sdf"
ACTIVES_FINAL_GZ = ALDR_DIR / "actives_final.sdf.gz"
ACTIVES_EXTENDED20_SDF = ALDR_DIR / "actives_extended20.sdf"
PROTEIN_LIST_JSON = AFDB_MOUSE_DIR / "protein_list.json"

METRICS_JSONL = OUTPUT_DIR / "nxm_poc_parallel_metrics.jsonl"
SUMMARY_TXT = OUTPUT_DIR / "nxm_poc_parallel_summary.txt"
HDF5_DIR = OUTPUT_DIR / "nxm_poc_hdf5"  # nxm_poc.py と同じリポジトリ（冪等性共有）

# pybel先にimportしてからvina後という順序を維持するため
# nxm_poc.py からユーティリティ関数を import
from nxm_poc import grid_box_centroid, load_proteins, load_compounds, prepare_extended20


def dock_one_protein(
    protein_pdb_path: str,
    protein_id: str,
    compound_sdf_path: str,
    compound_indices: List[int],
    grid_center: List[float],
    grid_size: List[float],
    exhaustiveness: int = 1,
    dry_run: bool = False,
) -> List[dict]:
    """1タンパク質の指定化合物群ドッキングを実行する。

    HDF5は一切開かない。結果をdictのリストとして返す。
    """
    import os
    import tempfile
    import time
    from pathlib import Path

    os.environ["BABEL_QUIET"] = "1"

    if not compound_indices:
        return []

    if dry_run:
        return [
            {
                "score": None,
                "sdf_content": None,
                "protein_content_hash": None,
                "compound_content_hash": None,
                "compoundset_content_hash": None,
                "compound_index": idx,
                "protein_id": protein_id,
                "compound_set_id": None,
                "metadata": {},
                "reused": False,
                "error": None,
                "dock_sec": 0.0,
            }
            for idx in compound_indices
        ]

    # pybel先にimport（順序厳守）
    from docking_automation.converters.molecule_converter import MoleculeConverter
    from docking_automation.docking.autodockvina_docking import AutoDockVina
    from docking_automation.molecule.compound_set import CompoundSet
    from docking_automation.molecule.protein import Protein
    from vina import Vina

    docking_tool = AutoDockVina()
    converter = MoleculeConverter()

    # タンパク質前処理
    protein = Protein(Path(protein_pdb_path), id=protein_id)
    try:
        prep_protein = docking_tool._preprocess_protein(protein)
    except Exception as e:
        return [
            {
                "score": None,
                "sdf_content": None,
                "protein_content_hash": protein.content_hash,
                "compound_content_hash": None,
                "compoundset_content_hash": None,
                "compound_index": idx,
                "protein_id": protein_id,
                "compound_set_id": None,
                "metadata": {},
                "reused": False,
                "error": f"タンパク質前処理失敗: {e}",
                "dock_sec": 0.0,
            }
            for idx in compound_indices
        ]

    protein_content_hash = protein.content_hash

    # 化合物セット読み込み・前処理（with_indicesはフィルタ不能のため全件前処理）
    # compound_indicesはoriginal_idxとしてprep_compounds.file_paths[original_idx]で参照
    full_compound_set = CompoundSet(Path(compound_sdf_path))
    try:
        prep_compounds = docking_tool._preprocess_compound_set(full_compound_set)
    except Exception as e:
        return [
            {
                "score": None,
                "sdf_content": None,
                "protein_content_hash": protein_content_hash,
                "compound_content_hash": None,
                "compoundset_content_hash": None,
                "compound_index": idx,
                "protein_id": protein_id,
                "compound_set_id": None,
                "metadata": {},
                "reused": False,
                "error": f"化合物前処理失敗: {e}",
                "dock_sec": 0.0,
            }
            for idx in compound_indices
        ]

    compoundset_content_hash = prep_compounds.content_hash

    # Vinaインスタンスをタンパク質毎に1回生成
    v = Vina(cpu=1, seed=1, verbosity=0)
    v.set_receptor(str(prep_protein.file_path))
    # compute_vina_maps もタンパク質毎に1回
    v.compute_vina_maps(
        center=[float(c) for c in grid_center],
        box_size=[float(s) for s in grid_size],
    )

    results = []
    for original_idx in compound_indices:
        t_dock_start = time.time()
        try:
            # with_indices が実装上フィルタ不能のため original_idx で直接参照
            compound_path = prep_compounds.file_paths[original_idx]
            compound_hash = prep_compounds.get_compound_hash(original_idx)

            temp_dir = Path(tempfile.mkdtemp())
            output_pdbqt = temp_dir / f"output_{original_idx}.pdbqt"
            output_sdf = temp_dir / f"output_{original_idx}.sdf"

            # 化合物ループ内: set_ligand → dock → energies のみ
            v.set_ligand_from_file(str(compound_path))
            v.dock(exhaustiveness=exhaustiveness, n_poses=3, min_rmsd=1.0)
            v.write_poses(str(output_pdbqt), n_poses=3, overwrite=True)
            converter.pdbqt_to_sdf(output_pdbqt, output_sdf)
            scores = v.energies()
            dock_sec = time.time() - t_dock_start

            results.append(
                {
                    "score": float(scores[0, 0]),
                    "sdf_content": output_sdf.read_text(),
                    "protein_content_hash": protein_content_hash,
                    "compound_content_hash": compound_hash,
                    "compoundset_content_hash": compoundset_content_hash,
                    "compound_index": original_idx,
                    "protein_id": protein_id,
                    "compound_set_id": compound_path.stem.split("_")[0],
                    "metadata": {
                        "tool": "AutoDock Vina",
                        "exhaustiveness": exhaustiveness,
                        "num_modes": 3,
                    },
                    "reused": False,
                    "error": None,
                    "dock_sec": round(dock_sec, 3),
                }
            )
        except Exception as e:
            dock_sec = time.time() - t_dock_start
            results.append(
                {
                    "score": None,
                    "sdf_content": None,
                    "protein_content_hash": protein_content_hash,
                    "compound_content_hash": None,
                    "compoundset_content_hash": compoundset_content_hash,
                    "compound_index": original_idx,
                    "protein_id": protein_id,
                    "compound_set_id": None,
                    "metadata": {},
                    "reused": False,
                    "error": str(e),
                    "dock_sec": round(dock_sec, 3),
                }
            )

    return results


def _save_result_dict(result_dict: dict, repo) -> None:
    """結果dictをDockingResultに変換してHDF5に保存する（メインプロセスのみ）。"""
    from docking_automation.docking.docking_result import DockingResult

    temp_dir = Path(tempfile.mkdtemp())
    output_sdf = temp_dir / f"result_{result_dict['compound_index']}.sdf"
    output_sdf.write_text(result_dict["sdf_content"])

    result = DockingResult(
        result_path=output_sdf,
        protein_id=result_dict["protein_id"],
        compound_set_id=result_dict["compound_set_id"],
        compound_index=result_dict["compound_index"],
        docking_score=result_dict["score"],
        protein_content_hash=result_dict["protein_content_hash"],
        compound_content_hash=result_dict["compound_content_hash"],
        compoundset_content_hash=result_dict["compoundset_content_hash"],
        metadata=result_dict["metadata"],
    )
    repo.save(result)


def run_nxm_parallel(
    proteins: List,
    compound_set,
    repo,
    metrics_fp,
    compound_sdf_path: Path,
    n_workers: int = 4,
    exhaustiveness: int = 1,
    dry_run: bool = False,
) -> List[dict]:
    """Dask LocalCluster を使った並列N×Mドッキング。"""
    from dask.distributed import Client, LocalCluster, as_completed

    from docking_automation.docking.autodockvina_docking import AutoDockVina

    docking_tool = AutoDockVina()

    # 化合物セットをメインプロセスで前処理してハッシュを取得
    prep_compounds = docking_tool._preprocess_compound_set(compound_set)
    n_compounds = len(prep_compounds.file_paths)
    all_compound_hashes = [prep_compounds.get_compound_hash(i) for i in range(n_compounds)]

    all_metrics: List[dict] = []

    with LocalCluster(
        n_workers=n_workers,
        threads_per_worker=1,
        processes=True,
        memory_limit="4GB",
    ) as cluster, Client(cluster) as client:

        futures_to_protein: dict = {}

        for protein in proteins:
            protein_id = protein.id

            # グリッドボックス計算
            try:
                grid_box = grid_box_centroid(protein.path)
            except Exception as e:
                print(f"[ERROR] {protein_id}: グリッドボックス計算失敗: {e}")
                continue

            center = [float(c) for c in grid_box.center]
            size = [float(s) for s in grid_box.size]

            # 未計算ペアをフィルタ（メインプロセスで事前チェック）
            unprocessed_indices: List[int] = []
            for idx in range(n_compounds):
                compound_hash = all_compound_hashes[idx]
                if repo._exists(protein.content_hash, compound_hash):
                    # 再利用済み：メトリクスに記録してスキップ
                    try:
                        result = repo.load_by_hashes(protein.content_hash, compound_hash)
                        score = result.docking_score if result is not None else None
                    except Exception:
                        score = None
                    metric = {
                        "protein_id": protein_id,
                        "compound_idx": idx,
                        "score": score,
                        "dock_sec": 0.0,
                        "reused": True,
                        "error": None,
                        "timestamp": int(time.time()),
                    }
                    metrics_fp.write(json.dumps(metric, ensure_ascii=False) + "\n")
                    metrics_fp.flush()
                    all_metrics.append(metric)
                    print(f"  [{protein_id}] compound[{idx}]: reused")
                else:
                    unprocessed_indices.append(idx)

            if not unprocessed_indices:
                print(f"[{protein_id}] 全ペア再利用済み。スキップ。")
                continue

            print(f"[{protein_id}] {len(unprocessed_indices)}件をキューに追加...")
            future = client.submit(
                dock_one_protein,
                str(protein.path),
                protein_id,
                str(compound_sdf_path),
                unprocessed_indices,
                center,
                size,
                exhaustiveness,
                dry_run,
                retries=2,
            )
            futures_to_protein[future] = protein_id

        # 結果収集・HDF5保存（メインプロセスで逐次）
        for future in as_completed(list(futures_to_protein.keys())):
            protein_id = futures_to_protein[future]
            try:
                results = future.result()
                for result_dict in results:
                    if result_dict.get("error"):
                        print(
                            f"  [ERROR] {protein_id} compound[{result_dict['compound_index']}]:"
                            f" {result_dict['error']}"
                        )
                        metric = {
                            "protein_id": protein_id,
                            "compound_idx": result_dict["compound_index"],
                            "score": None,
                            "dock_sec": result_dict.get("dock_sec", 0.0),
                            "reused": False,
                            "error": result_dict["error"],
                            "timestamp": int(time.time()),
                        }
                    else:
                        if not dry_run and result_dict.get("sdf_content"):
                            _save_result_dict(result_dict, repo)
                        score = result_dict["score"]
                        dock_sec = result_dict.get("dock_sec", 0.0)
                        status = "dry-run" if dry_run else f"score={score:.2f}"
                        print(
                            f"  [{protein_id}] compound[{result_dict['compound_index']}]:"
                            f" {status} ({dock_sec:.1f}s)"
                        )
                        metric = {
                            "protein_id": protein_id,
                            "compound_idx": result_dict["compound_index"],
                            "score": score,
                            "dock_sec": dock_sec,
                            "reused": False,
                            "error": None,
                            "timestamp": int(time.time()),
                        }

                    metrics_fp.write(json.dumps(metric, ensure_ascii=False) + "\n")
                    metrics_fp.flush()
                    all_metrics.append(metric)

            except Exception as e:
                print(f"[ERROR] {protein_id} タスク失敗: {e}")

    return all_metrics


def write_summary_parallel(all_metrics: List[dict], elapsed_total: float) -> None:
    """並列版サマリファイルを書き出す。"""
    if not all_metrics:
        with open(SUMMARY_TXT, "w") as f:
            f.write("メトリクスなし（エラーで全ペア失敗）\n")
        return

    n_total = len(all_metrics)
    n_reused = sum(1 for m in all_metrics if m["reused"])
    n_new = n_total - n_reused
    scores = [m["score"] for m in all_metrics if m["score"] is not None]
    avg_score = np.mean(scores) if scores else None
    new_metrics = [m for m in all_metrics if not m["reused"]]
    avg_dock = np.mean([m["dock_sec"] for m in new_metrics]) if new_metrics else 0.0

    lines = [
        "=== N×M ドッキングPoC (Dask並列版) サマリ ===",
        f"総ペア数:        {n_total}",
        f"新規ドッキング:  {n_new}",
        f"再利用:          {n_reused}",
        f"総経過時間:      {elapsed_total:.1f}s",
        (
            f"平均ドッキング時間: {avg_dock:.2f}s/ペア (新規のみ)"
            if n_new > 0
            else "平均ドッキング時間: N/A"
        ),
        (
            f"平均スコア:      {avg_score:.3f}"
            if avg_score is not None
            else "平均スコア: N/A"
        ),
    ]
    text = "\n".join(lines) + "\n"
    print("\n" + text)
    with open(SUMMARY_TXT, "w") as f:
        f.write(text)


def parse_args():
    parser = argparse.ArgumentParser(description="N×M ドッキングPoC (Dask並列版)")
    parser.add_argument(
        "--workers",
        type=int,
        default=min(os.cpu_count() or 4, 4),
        help="Daskワーカー数（デフォルト: min(cpu_count, 4)）",
    )
    parser.add_argument(
        "--exhaustiveness",
        type=int,
        default=1,
        help="Vina exhaustiveness（デフォルト: 1）",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="設定確認のみ（ドッキング実行なし）",
    )
    parser.add_argument(
        "--max-compounds",
        type=int,
        default=10,
        help="処理する化合物の最大数（デフォルト: 10）",
    )
    parser.add_argument(
        "--extended",
        action="store_true",
        help="actives_extended20.sdf を使用（20件）",
    )
    parser.add_argument(
        "--protein-list",
        type=Path,
        default=None,
        help="protein_list.json のパス（デフォルト: afdb_mouse/protein_list.json）",
    )
    return parser.parse_args()


def main():
    args = parse_args()

    print("=== N×M ドッキングPoC (Dask並列版) 開始 ===")
    print(f"ワーカー数: {args.workers}")
    print(f"作業ディレクトリ: {PROJECT_DIR}")

    # actives_extended20.sdf 準備
    if ACTIVES_FINAL_GZ.exists():
        prepare_extended20()

    # 入力ファイル選択
    if args.extended:
        compound_sdf = ACTIVES_EXTENDED20_SDF
        max_n = args.max_compounds if args.max_compounds != 10 else 20
    else:
        compound_sdf = ACTIVES_SUBSET_SDF
        max_n = args.max_compounds

    if not compound_sdf.exists():
        print(f"[ERROR] 化合物ファイルが見つかりません: {compound_sdf}")
        sys.exit(1)

    proteins = load_proteins(args.protein_list)
    if not proteins:
        print("[ERROR] タンパク質が1件も読み込めませんでした。")
        sys.exit(1)

    compound_set = load_compounds(compound_sdf, max_n=max_n)

    print(f"\n実行設定:")
    print(f"  タンパク質数:   {len(proteins)}")
    print(f"  化合物数:       {max_n}")
    print(f"  総ペア数:       {len(proteins) * max_n}")
    print(f"  exhaustiveness: {args.exhaustiveness}")
    print(f"  Daskワーカー数: {args.workers}")
    print(f"  HDF5 出力先:    {HDF5_DIR}")

    if args.dry_run:
        print("\n[dry-run] 設定確認完了。ドッキングは実行しません。")
        return

    # 出力ディレクトリ準備
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    HDF5_DIR.mkdir(parents=True, exist_ok=True)

    # HDF5 リポジトリ初期化
    from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
        DockingResultRepositoryFactory,
        RepositoryType,
    )

    repo = DockingResultRepositoryFactory.create(
        RepositoryType.HDF5,
        base_directory=HDF5_DIR,
        config={"mode": "append"},
    )

    # N×M ドッキング実行（Dask並列）
    t_start = time.time()
    with open(METRICS_JSONL, "w") as metrics_fp:
        all_metrics = run_nxm_parallel(
            proteins=proteins,
            compound_set=compound_set,
            repo=repo,
            metrics_fp=metrics_fp,
            compound_sdf_path=compound_sdf,
            n_workers=args.workers,
            exhaustiveness=args.exhaustiveness,
            dry_run=args.dry_run,
        )
    elapsed = time.time() - t_start

    # サマリ出力
    write_summary_parallel(all_metrics, elapsed)
    print(f"\nメトリクス: {METRICS_JSONL}")
    print(f"サマリ:     {SUMMARY_TXT}")
    print("=== 完了 ===")


if __name__ == "__main__":
    main()
