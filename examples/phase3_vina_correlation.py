"""Phase 3 Wave 5: Vina vs UniDock score correlation (10x10 subset).

設計書 Q3: Pearson r > 0.9 を確認する。
新規ディレクトリを使用してHDF5キャッシュを避ける。

実行方法:
  cd /workspaces/20260422_mouse_docking/docking_automation
  python examples/phase3_vina_correlation.py 2>&1 | tee /tmp/phase3_corr.log
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr

project_root = Path(__file__).parent.parent
sys.path.insert(0, str(project_root))

from docking_automation.docking.fpocket_grid_box_predictor import FpocketGridBoxPredictor
from docking_automation.docking.grid_box_cache import GridBoxCache
from docking_automation.docking.screening_runner import ScreeningRunner
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein
from docking_automation.molecule.protein_set import ProteinSet

DATA_DIR = project_root / "examples" / "input"
PROTEIN_LIST_PATH = DATA_DIR / "afdb_mouse" / "protein_list.json"
COMPOUND_SDF_PATH = DATA_DIR / "ALDR" / "actives_subset.sdf"

OUTPUT_DIR = project_root / "examples" / "output"
CORR_DIR = Path("/tmp/corr_wave5")
VINA_HDF5 = CORR_DIR / "vina.h5"
UNIDOCK_HDF5 = CORR_DIR / "unidock.h5"
GRID_BOX_CACHE_PATH = CORR_DIR / "grid_cache.json"

N_PROTEINS = 10


def build_protein_set() -> ProteinSet:
    with open(PROTEIN_LIST_PATH) as f:
        protein_paths = json.load(f)
    protein_list = [Protein(Path(p)) for p in protein_paths[:N_PROTEINS]]
    return ProteinSet(protein_list)


def build_grid_box_cache(protein_set: ProteinSet) -> GridBoxCache:
    if GRID_BOX_CACHE_PATH.exists():
        grid_box_cache = GridBoxCache.from_file(GRID_BOX_CACHE_PATH)
        print("  GridBoxキャッシュ読み込み済み")
        return grid_box_cache
    grid_box_cache = GridBoxCache(GRID_BOX_CACHE_PATH)

    predictor = FpocketGridBoxPredictor()
    built = 0
    for protein in protein_set:
        try:
            grid_box = predictor.predict(protein, pocket_rank=1)
            grid_box_cache.put(protein, grid_box, source="fpocket")
            print(f"  [OK] {protein.id}")
            built += 1
        except Exception as e:
            print(f"  [NG] {protein.id}: {e}")
    grid_box_cache.save(GRID_BOX_CACHE_PATH)
    print(f"  GridBox構築完了: {built}/{N_PROTEINS}")
    return grid_box_cache


def read_scores_from_hdf5(hdf5_path: Path) -> dict[tuple[str, str], float]:
    """HDF5から (protein_hash, compound_hash) -> score の辞書を返す。"""
    import h5py
    scores: dict[tuple[str, str], float] = {}
    with h5py.File(hdf5_path, "r") as f:
        results_group = f.get("results", {})
        for p_hash in results_group:
            for c_hash in results_group[p_hash]:
                ds = results_group[p_hash][c_hash].get("docking_score")
                if ds is not None:
                    val = float(ds[()])
                    scores[(p_hash, c_hash)] = val
    return scores


def run_backend(
    backend: str,
    protein_set: ProteinSet,
    compound_set: CompoundSet,
    grid_box_cache: GridBoxCache,
    hdf5_path: Path,
    search_mode: str = "balance",
) -> ScreeningRunner:
    kwargs: dict = dict(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=grid_box_cache,
        hdf5_path=hdf5_path,
        dask_n_workers=4,
        exhaustiveness=1,
        log_path=hdf5_path.with_suffix(".jsonl"),
        grid_box_missing_policy="skip",
        backend=backend,
    )
    if backend == "unidock":
        kwargs["search_mode"] = search_mode

    runner = ScreeningRunner(**kwargs)
    t0 = time.time()
    result = runner.run(resume=False)
    elapsed = time.time() - t0
    print(f"  {backend}: total={result.total_pairs}, new={result.new_pairs}, "
          f"failed={result.failed_pairs}, elapsed={elapsed:.1f}s")
    return runner, result


def main() -> None:
    CORR_DIR.mkdir(parents=True, exist_ok=True)

    # ===== 1. データロード =====
    print("=== Step 1: タンパク質・化合物ロード ===")
    protein_set = build_protein_set()
    compound_set = CompoundSet.create(COMPOUND_SDF_PATH)
    n_compounds = compound_set.get_compound_count()
    print(f"  タンパク質: {len(protein_set)} 件, 化合物: {n_compounds} 件")

    # ===== 2. GridBox構築 =====
    print("\n=== Step 2: GridBoxCache 構築 ===")
    grid_box_cache = build_grid_box_cache(protein_set)

    # ===== 3. Vina 実行 =====
    print("\n=== Step 3: Vina 実行 ===")
    _, vina_result = run_backend(
        backend="vina",
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=grid_box_cache,
        hdf5_path=VINA_HDF5,
    )

    # ===== 4. UniDock 実行 =====
    print("\n=== Step 4: UniDock 実行 (search_mode=balance) ===")
    _, unidock_result = run_backend(
        backend="unidock",
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=grid_box_cache,
        hdf5_path=UNIDOCK_HDF5,
        search_mode="balance",
    )

    # ===== 5. スコア抽出と相関計算 =====
    print("\n=== Step 5: スコア相関計算 ===")
    vina_scores_map = read_scores_from_hdf5(VINA_HDF5)
    unidock_scores_map = read_scores_from_hdf5(UNIDOCK_HDF5)

    common_keys = set(vina_scores_map.keys()) & set(unidock_scores_map.keys())
    print(f"  Vina スコアペア数: {len(vina_scores_map)}")
    print(f"  UniDock スコアペア数: {len(unidock_scores_map)}")
    print(f"  共通ペア数: {len(common_keys)}")

    # NOTE: dock_one_protein() はペナルティフィルタを未適用のため生スコアに正値が混入。
    # DockingParameters デフォルト閾値 (-30 < score < 5) で両 backend を揃えてフィルタする。
    SCORE_MIN, SCORE_MAX = -30.0, 5.0
    valid_keys = [
        k for k in common_keys
        if SCORE_MIN < vina_scores_map[k] < SCORE_MAX
        and SCORE_MIN < unidock_scores_map[k] < SCORE_MAX
    ]
    raw_total = len(common_keys)
    print(f"  ペナルティ除外後 有効ペア: {len(valid_keys)}/{raw_total}")

    if len(valid_keys) < 2:
        print("[FAIL] 有効ペア数が不足しています (< 2)")
        sys.exit(1)

    vina_arr = np.array([vina_scores_map[k] for k in valid_keys])
    unidock_arr = np.array([unidock_scores_map[k] for k in valid_keys])

    r, p_val = pearsonr(vina_arr, unidock_arr)
    print(f"\n  Valid pairs:  {len(vina_arr)}/{raw_total}")
    print(f"  Pearson r =   {r:.4f} (p = {p_val:.4e})")
    print(f"  Vina scores:    mean={np.mean(vina_arr):.2f}, "
          f"range=[{vina_arr.min():.2f}, {vina_arr.max():.2f}]")
    print(f"  UniDock scores: mean={np.mean(unidock_arr):.2f}, "
          f"range=[{unidock_arr.min():.2f}, {unidock_arr.max():.2f}]")

    # ===== 6. Q3 判定 =====
    print("\n=== Step 6: Q3 判定 ===")
    q3_passed = r > 0.9
    if q3_passed:
        print(f"  [PASS] Q3: Pearson r = {r:.4f} > 0.9 達成")
    else:
        print(f"  [WARN] Q3 未達: r = {r:.4f} < 0.9")
        print("  search_mode='detail' での再実行を検討してください")

    print(f"\n  [BUG報告] dock_one_protein() にペナルティフィルタが未適用。")
    print(f"  正値スコア混入: Vina {sum(1 for k in common_keys if vina_scores_map[k] >= SCORE_MAX)}件, "
          f"UniDock {sum(1 for k in common_keys if unidock_scores_map[k] >= SCORE_MAX)}件 / {raw_total}件")
    print(f"  修正案: screening_runner.py dock_one_protein() の UniDock/Vina 両パスに閾値フィルタ追加が必要。")

    assert q3_passed, f"Pearson r = {r:.4f} < 0.9 (設計書 Q3 未達)"

    return {
        "valid_pairs": len(vina_arr),
        "pearson_r": round(float(r), 4),
        "q3_passed": q3_passed,
        "vina_mean": round(float(np.mean(vina_arr)), 2),
        "unidock_filtered_mean": round(float(np.mean(unidock_arr)), 2),
    }


if __name__ == "__main__":
    main()
