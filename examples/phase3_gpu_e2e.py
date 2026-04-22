"""Phase 3 GPU E2E: ScreeningRunner(backend='unidock') で 10x100 GPU ドッキング."""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

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
COMPOUND_SDF_PATH = DATA_DIR / "ALDR" / "actives_subset100.sdf"
OUTPUT_DIR = project_root / "examples" / "output"
HDF5_PATH = OUTPUT_DIR / "phase3_gpu_e2e.h5"
GRID_BOX_CACHE_PATH = OUTPUT_DIR / "phase3_gpu_e2e_grid_cache.json"
LOG_PATH = OUTPUT_DIR / "phase3_gpu_e2e_screening.jsonl"

N_PROTEINS = 10
N_COMPOUNDS = 100


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ===== 1. タンパク質 10 件をロード =====
    print("=== Step 1: タンパク質 10 件をロード ===")
    with open(PROTEIN_LIST_PATH) as f:
        protein_paths = json.load(f)

    protein_list = [Protein(Path(p)) for p in protein_paths[:N_PROTEINS]]
    protein_set = ProteinSet(protein_list)
    print(f"  タンパク質数: {len(protein_set)}")
    for p in protein_set:
        print(f"  - {p.id}")

    # ===== 2. 化合物 100 件をロード =====
    print(f"\n=== Step 2: 化合物 {N_COMPOUNDS} 件をロード ===")
    compound_set = CompoundSet.create(COMPOUND_SDF_PATH)
    compound_count = compound_set.get_compound_count()
    print(f"  化合物数: {compound_count}")
    assert compound_count == N_COMPOUNDS, f"化合物数不一致: {compound_count} != {N_COMPOUNDS}"

    # ===== 3. GridBoxCache 構築 (fpocket) =====
    print("\n=== Step 3: GridBoxCache 構築 (fpocket) ===")
    grid_box_cache = GridBoxCache(GRID_BOX_CACHE_PATH)
    predictor = FpocketGridBoxPredictor()

    built_count = 0
    skipped_count = 0
    for protein in protein_set:
        t0 = time.monotonic()
        try:
            grid_box = predictor.predict(protein, pocket_rank=1)
            grid_box_cache.put(protein, grid_box, source="fpocket")
            elapsed = time.monotonic() - t0
            print(f"  [OK] {protein.id}: center={grid_box.center.tolist()}, elapsed={elapsed:.2f}s")
            built_count += 1
        except Exception as e:
            elapsed = time.monotonic() - t0
            print(f"  [NG] {protein.id}: {e} (elapsed={elapsed:.2f}s)")
            skipped_count += 1

    print(f"  GridBox 構築完了: {built_count} 件 / スキップ: {skipped_count} 件")
    grid_box_cache.save(GRID_BOX_CACHE_PATH)

    # ===== 4. ScreeningRunner 作成 (backend=unidock) =====
    print("\n=== Step 4: ScreeningRunner 初期化 (backend=unidock) ===")
    runner = ScreeningRunner(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=grid_box_cache,
        hdf5_path=HDF5_PATH,
        dask_n_workers=4,
        exhaustiveness=1,
        log_path=LOG_PATH,
        grid_box_missing_policy="skip",
        backend="unidock",
        search_mode="balance",
    )
    print(f"  backend: {runner.backend}")
    print(f"  search_mode: {runner.search_mode}")
    print(f"  total_pairs予定: {len(protein_set) * compound_count}")

    # ===== 5. Run 1 =====
    print("\n=== Step 5: Run 1 (初回 GPU ドッキング) ===")
    t_run1 = time.time()
    result1 = runner.run(resume=True)
    elapsed1 = time.time() - t_run1

    print("  Run1 結果:")
    print(f"    total_pairs:  {result1.total_pairs}")
    print(f"    new_pairs:    {result1.new_pairs}")
    print(f"    reused_pairs: {result1.reused_pairs}")
    print(f"    failed_pairs: {result1.failed_pairs}")
    print(f"    elapsed:      {elapsed1:.1f}s")

    assert result1.failed_pairs == 0, (
        f"GPU E2E failed: {result1.failed_pairs} ペアが失敗"
    )

    # ===== 6. Run 2 (冪等性検証) =====
    print("\n=== Step 6: Run 2 (冪等性検証) ===")
    t_run2 = time.time()
    result2 = runner.run(resume=True)
    elapsed2 = time.time() - t_run2

    print("  Run2 結果:")
    print(f"    total_pairs:  {result2.total_pairs}")
    print(f"    new_pairs:    {result2.new_pairs}")
    print(f"    reused_pairs: {result2.reused_pairs}")
    print(f"    failed_pairs: {result2.failed_pairs}")
    print(f"    elapsed:      {elapsed2:.2f}s")

    assert result2.new_pairs == 0, (
        f"冪等性違反: Run2.new_pairs == {result2.new_pairs} (期待値: 0)"
    )
    assert result2.reused_pairs == result1.new_pairs, (
        f"Run2.reused_pairs == {result2.reused_pairs} (期待値: {result1.new_pairs})"
    )

    # ===== 7. スコア統計 =====
    print("\n=== Step 7: スコア統計 ===")
    import h5py
    scores = []
    if HDF5_PATH.exists():
        with h5py.File(HDF5_PATH, "r") as f:
            results_group = f.get("results", {})
            for p_hash in results_group:
                for c_hash in results_group[p_hash]:
                    ds = results_group[p_hash][c_hash].get("docking_score")
                    if ds is not None:
                        scores.append(float(ds[()]))

    if scores:
        print(f"  スコア件数: {len(scores)}")
        print(f"  Score range: {min(scores):.2f} ~ {max(scores):.2f} kcal/mol")
        print(f"  Mean score:  {sum(scores)/len(scores):.2f} kcal/mol")
    else:
        print("  スコアデータなし")

    # ===== 8. Phase 2 比較 =====
    print("\n=== Step 8: Phase 2 比較 ===")
    phase2_hdf5 = OUTPUT_DIR / "phase2_e2e.h5"
    if phase2_hdf5.exists():
        size2 = phase2_hdf5.stat().st_size
        size3 = HDF5_PATH.stat().st_size if HDF5_PATH.exists() else 0
        print(f"  Phase 2 HDF5: {size2/1024:.1f} KB ({size2} bytes, 10×10=100 ペア)")
        print(f"  Phase 3 HDF5: {size3/1024:.1f} KB ({size3} bytes, 10×{N_COMPOUNDS}={result1.new_pairs} ペア)")
    else:
        print("  Phase 2 HDF5 なし (比較不可)")

    # ===== 9. 完了サマリー =====
    print("\n=== 完了 ===")
    print(f"  Run1: new={result1.new_pairs}, reused={result1.reused_pairs}, failed={result1.failed_pairs}, elapsed={elapsed1:.1f}s")
    print(f"  Run2: new={result2.new_pairs}, reused={result2.reused_pairs}, failed={result2.failed_pairs}, elapsed={elapsed2:.2f}s")
    if scores:
        print(f"  Score: mean={sum(scores)/len(scores):.2f}, min={min(scores):.2f}, max={max(scores):.2f} kcal/mol")
    print("  [PASS] GPU E2E 完了")

    return result1, result2, scores


if __name__ == "__main__":
    main()
