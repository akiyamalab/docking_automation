"""Phase 2 E2E 動作確認: 10タンパク質 × 10化合物

実行方法:
  cd /workspaces/20260422_mouse_docking/docking_automation
  python examples/phase2_e2e.py
"""
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
COMPOUND_SDF_PATH = DATA_DIR / "ALDR" / "actives_subset.sdf"
OUTPUT_DIR = project_root / "examples" / "output"
HDF5_PATH = OUTPUT_DIR / "phase2_e2e.h5"
GRID_BOX_CACHE_PATH = OUTPUT_DIR / "phase2_e2e_grid_cache.json"
LOG_PATH = OUTPUT_DIR / "phase2_e2e_screening.jsonl"


def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    # ===== 1. タンパク質 10 件を ProteinSet に =====
    print("=== Step 1: タンパク質 10 件をロード ===")
    with open(PROTEIN_LIST_PATH) as f:
        protein_paths = json.load(f)

    protein_list = [Protein(Path(p)) for p in protein_paths[:10]]
    protein_set = ProteinSet(protein_list)
    print(f"  タンパク質数: {len(protein_set)}")
    for p in protein_set:
        print(f"  - {p.id}")

    # ===== 2. 化合物 10 件をロード =====
    print("\n=== Step 2: 化合物 10 件をロード ===")
    compound_set = CompoundSet.create(COMPOUND_SDF_PATH)
    compound_count = compound_set.get_compound_count()
    print(f"  化合物数: {compound_count}")

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

    # ===== 4. ScreeningRunner 作成 =====
    print("\n=== Step 4: ScreeningRunner 初期化 ===")
    runner = ScreeningRunner(
        protein_set=protein_set,
        compound_set=compound_set,
        grid_box_cache=grid_box_cache,
        hdf5_path=HDF5_PATH,
        dask_n_workers=4,
        exhaustiveness=1,
        log_path=LOG_PATH,
        grid_box_missing_policy="skip",
    )
    print(f"  total_pairs予定: {len(protein_set) * compound_count}")

    # ===== 5. Run 1 =====
    print("\n=== Step 5: Run 1 (初回) ===")
    result1 = runner.run(resume=True)

    print("  Run1 結果:")
    print(f"    total_pairs:  {result1.total_pairs}")
    print(f"    new_pairs:    {result1.new_pairs}")
    print(f"    reused_pairs: {result1.reused_pairs}")
    print(f"    failed_pairs: {result1.failed_pairs}")
    print(f"    elapsed:      {result1.elapsed_sec:.2f}s")

    # ===== 6. Run 2 (resume 検証) =====
    print("\n=== Step 6: Run 2 (resume検証) ===")
    result2 = runner.run(resume=True)

    print("  Run2 結果:")
    print(f"    total_pairs:  {result2.total_pairs}")
    print(f"    new_pairs:    {result2.new_pairs}")
    print(f"    reused_pairs: {result2.reused_pairs}")
    print(f"    failed_pairs: {result2.failed_pairs}")
    print(f"    elapsed:      {result2.elapsed_sec:.2f}s")

    # ===== 7. 検証 =====
    print("\n=== Step 7: 検証 ===")
    errors = []

    if result1.new_pairs == 0 and built_count > 0:
        errors.append(f"Run1.new_pairs == 0 だが GridBox は {built_count} 件構築されている")

    if result2.new_pairs != 0:
        errors.append(f"Run2.new_pairs == {result2.new_pairs} (期待値: 0 = 冪等性違反)")

    expected_reused = result1.new_pairs
    if result2.reused_pairs != expected_reused:
        errors.append(
            f"Run2.reused_pairs == {result2.reused_pairs} (期待値: {expected_reused})"
        )

    if not HDF5_PATH.exists():
        errors.append(f"HDF5ファイルが存在しない: {HDF5_PATH}")

    if errors:
        print("  [FAIL] 検証エラー:")
        for e in errors:
            print(f"    - {e}")
    else:
        print("  [PASS] 全検証クリア!")
        print(f"  [PASS] Run2.new_pairs == 0 (冪等性確認)")
        print(f"  [PASS] Run2.reused_pairs == {result2.reused_pairs} (reuse確認)")

    # ===== 8. HDF5 サイズ計測 =====
    print("\n=== Step 8: HDF5 サイズ計測 ===")
    if HDF5_PATH.exists():
        size_bytes = HDF5_PATH.stat().st_size
        size_kb = size_bytes / 1024
        print(f"  HDF5サイズ: {size_kb:.1f} KB ({size_bytes} bytes)")

        actual_pairs = result1.new_pairs
        if actual_pairs > 0:
            bytes_per_pair = size_bytes / actual_pairs
            extrap_bytes = bytes_per_pair * 2e8
            extrap_gb = extrap_bytes / (1024 ** 3)
            print(f"  bytes/pair: {bytes_per_pair:.1f}")
            print(f"  Phase 4 外挿 (2e8 ペア): {extrap_gb:.1f} GB")
        else:
            print("  外挿不可 (new_pairs == 0)")
    else:
        print("  HDF5ファイルなし")

    print("\n=== 完了 ===")
    print(f"  Run1: new={result1.new_pairs}, reused={result1.reused_pairs}, failed={result1.failed_pairs}, elapsed={result1.elapsed_sec:.2f}s")
    print(f"  Run2: new={result2.new_pairs}, reused={result2.reused_pairs}, failed={result2.failed_pairs}, elapsed={result2.elapsed_sec:.2f}s")

    return result1, result2, built_count


if __name__ == "__main__":
    main()
