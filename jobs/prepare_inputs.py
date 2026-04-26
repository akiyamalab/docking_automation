"""N×M docking 入力前処理 (unidock2 向け).

N 受容体 (PDB) + M リガンド (SDF を 1 件ずつ分割) + N グリッドボックス を
jobs/inputs/ に書き出す。

unidock2 は PDB/SDF をそのまま受けるので PDBQT 変換は不要。

Usage (from /workspaces/20260422_mouse_docking/docking_automation):
    python jobs/prepare_inputs.py --n-receptors 10  --n-ligands 10
    python jobs/prepare_inputs.py --n-receptors 100 --n-ligands 100
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

DOCKING_AUTOMATION = Path("/workspaces/20260422_mouse_docking/docking_automation")
sys.path.insert(0, str(DOCKING_AUTOMATION))

EXAMPLES = DOCKING_AUTOMATION / "examples"
PROTEIN_LIST = EXAMPLES / "input" / "afdb_mouse" / "protein_list.json"
GRID_CACHE = EXAMPLES / "output" / "midscale_grid_cache.json"

LIGAND_SDF_10 = EXAMPLES / "input" / "ALDR" / "actives_subset.sdf"
LIGAND_SDF_100 = EXAMPLES / "input" / "ALDR" / "actives_subset100.sdf"

PROJECT_ROOT = Path("/workspaces/20260422_mouse_docking/docking_automation")
INPUTS_DIR = PROJECT_ROOT / "jobs" / "inputs"
RECEPTORS_DIR = INPUTS_DIR / "receptors_pdb"
LIGANDS_DIR = INPUTS_DIR / "ligands_sdf"
BOXES_JSON = INPUTS_DIR / "boxes.json"
BOXES_TSV = INPUTS_DIR / "boxes.tsv"


def _clear(d: Path, glob: str) -> None:
    if d.exists():
        for p in d.glob(glob):
            p.unlink()


def split_sdf(sdf_path: Path, out_dir: Path, n: int) -> int:
    """SDF を 1 分子ずつ分割。RDKit を使う。"""
    from rdkit import Chem
    suppl = Chem.SDMolSupplier(str(sdf_path))
    count = 0
    for i, mol in enumerate(suppl):
        if mol is None:
            continue
        if count >= n:
            break
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"lig_{i:04d}"
        # 安全なファイル名に
        safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in name)
        out = out_dir / f"{safe}.sdf"
        # 重複対策
        if out.exists():
            out = out_dir / f"{safe}_{i}.sdf"
        w = Chem.SDWriter(str(out))
        w.write(mol)
        w.close()
        count += 1
    return count


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--n-receptors", type=int, default=10)
    parser.add_argument("--n-ligands", type=int, default=10)
    args = parser.parse_args()
    n_rcp = args.n_receptors
    n_lig = args.n_ligands

    _clear(RECEPTORS_DIR, "*.pdb")
    _clear(LIGANDS_DIR, "*.sdf")
    RECEPTORS_DIR.mkdir(parents=True, exist_ok=True)
    LIGANDS_DIR.mkdir(parents=True, exist_ok=True)

    with open(PROTEIN_LIST) as f:
        protein_paths = [Path(p) for p in json.load(f)[:n_rcp]]
    if not protein_paths[0].is_absolute():
        protein_paths = [DOCKING_AUTOMATION / p for p in protein_paths]

    with open(GRID_CACHE) as f:
        grid_cache = json.load(f)["entries"]

    boxes = {}
    print(f"=== {len(protein_paths)} 受容体 PDB をコピー ===")
    skipped = []
    for i, pdb in enumerate(protein_paths, 1):
        rid = pdb.stem
        if rid not in grid_cache:
            skipped.append(rid)
            continue
        shutil.copy2(pdb, RECEPTORS_DIR / f"{rid}.pdb")
        if i % 10 == 0 or i == len(protein_paths):
            print(f"  [{i:3d}/{len(protein_paths)}] {rid}.pdb")
        e = grid_cache[rid]
        boxes[rid] = {"center": e["center"], "size": e["size"]}

    with open(BOXES_JSON, "w") as f:
        json.dump(boxes, f, indent=2)
    with open(BOXES_TSV, "w") as f:
        f.write("id\tcx\tcy\tcz\tsx\tsy\tsz\n")
        for rid, b in boxes.items():
            f.write(f"{rid}\t" + "\t".join(map(str, b["center"] + b["size"])) + "\n")
    print(f"  boxes: {len(boxes)} (skipped {len(skipped)})")

    sdf = LIGAND_SDF_100 if n_lig > 10 else LIGAND_SDF_10
    print(f"\n=== {n_lig} リガンド SDF 分割 (source: {sdf.name}) ===")
    written = split_sdf(sdf, LIGANDS_DIR, n_lig)
    print(f"  生成 SDF: {written} 件")

    print("\n=== 完了 ===")
    print(f"  受容体:  {RECEPTORS_DIR}/*.pdb ({len(boxes)} 件)")
    print(f"  リガンド: {LIGANDS_DIR}/*.sdf ({written} 件)")
    print(f"  ボックス: {BOXES_JSON} / {BOXES_TSV}")


if __name__ == "__main__":
    main()
