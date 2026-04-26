"""data/afdb/pdb/ から 1000 受容体を選び inputs/receptors_pdb/ にコピー.

DMS cache prep 用 (grid box 不要). 残基数 500 程度のフィルタを
ファイルサイズ < 400KB で代用 (PDB 1 残基 ≒ 800B のため).
"""
from __future__ import annotations

import argparse
import shutil
from pathlib import Path

DOCKING_AUTOMATION = Path("/workspaces/20260422_mouse_docking/docking_automation")
AFDB_TREE = DOCKING_AUTOMATION / "data" / "afdb" / "pdb"
INPUTS_DIR = DOCKING_AUTOMATION / "jobs" / "inputs"
RECEPTORS_DIR = INPUTS_DIR / "receptors_pdb"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n", type=int, default=1000)
    parser.add_argument("--max-size-kb", type=int, default=400,
                        help="残基数フィルタの代理 (1 残基 ≒ 800B)")
    args = parser.parse_args()

    if not AFDB_TREE.exists():
        raise SystemExit(f"AFDB tree not found: {AFDB_TREE}")

    # 既存の 100 件をクリア
    if RECEPTORS_DIR.exists():
        for p in RECEPTORS_DIR.glob("AF-*.pdb"):
            p.unlink()
    RECEPTORS_DIR.mkdir(parents=True, exist_ok=True)

    max_bytes = args.max_size_kb * 1024
    print(f"=== AFDB tree から最大 {args.n} 件選定 (size < {args.max_size_kb}KB) ===")

    selected: list[Path] = []
    # ソートして安定順序に
    for pdb in sorted(AFDB_TREE.rglob("AF-*.pdb")):
        if pdb.stat().st_size > max_bytes:
            continue
        selected.append(pdb)
        if len(selected) >= args.n:
            break

    print(f"  selected: {len(selected)}")
    for i, src in enumerate(selected, 1):
        shutil.copy2(src, RECEPTORS_DIR / src.name)
        if i % 100 == 0 or i == len(selected):
            print(f"  copied [{i:4d}/{len(selected)}]")

    print(f"=== 完了: {RECEPTORS_DIR} に {len(selected)} 件 ===")
    print(f"  total size: {sum(p.stat().st_size for p in RECEPTORS_DIR.glob('AF-*.pdb')) // 1024 // 1024} MB")


if __name__ == "__main__":
    main()
