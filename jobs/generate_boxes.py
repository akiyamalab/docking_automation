"""inputs/receptors_pdb/ にある全 PDB に対して fpocket で box を計算し
boxes.tsv / boxes.json を生成する.

- 既存の midscale_grid_cache.json に登録のある受容体は cache を使用
- それ以外は fpocket で予測 (rank 1)
- fpocket でも box が得られなければ centroid + 30Å にフォールバック

ローカル実行 (要 fpocket バイナリ + docking_automation Python 環境).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

DOCKING_AUTOMATION = Path("/workspaces/20260422_mouse_docking/docking_automation")
sys.path.insert(0, str(DOCKING_AUTOMATION))

from docking_automation.docking.fpocket_grid_box_predictor import FpocketGridBoxPredictor
from docking_automation.molecule.protein import Protein

GRID_CACHE = DOCKING_AUTOMATION / "examples" / "output" / "midscale_grid_cache.json"

INPUTS = Path("/workspaces/20260422_mouse_docking/docking_automation/jobs/inputs")
RECEPTORS = INPUTS / "receptors_pdb"
BOXES_JSON = INPUTS / "boxes.json"
BOXES_TSV = INPUTS / "boxes.tsv"
CACHE_LOCAL = INPUTS / "_fpocket_local_cache.json"  # 大規模なので increment cache

CENTROID_SIZE = 30.0


def ca_centroid(pdb: Path) -> tuple[float, float, float] | None:
    xs, ys, zs = [], [], []
    with open(pdb) as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                xs.append(float(line[30:38]))
                ys.append(float(line[38:46]))
                zs.append(float(line[46:54]))
            except ValueError:
                continue
    if not xs:
        return None
    return (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))


def predict_fpocket(pdb: Path, predictor: FpocketGridBoxPredictor) -> dict | None:
    try:
        protein = Protein(pdb)
        gb = predictor.predict(protein, pocket_rank=1)
        return {"center": gb.center.tolist(), "size": gb.size.tolist(), "source": "fpocket"}
    except Exception as e:
        print(f"  [warn] fpocket failed for {pdb.stem}: {type(e).__name__}: {e}")
        return None


def main() -> None:
    grid_cache = json.loads(GRID_CACHE.read_text())["entries"] if GRID_CACHE.exists() else {}
    local_cache = json.loads(CACHE_LOCAL.read_text()) if CACHE_LOCAL.exists() else {}

    pdbs = sorted(RECEPTORS.glob("AF-*.pdb"))
    print(f"=== {len(pdbs)} 受容体に対する box 生成 ===")
    print(f"  midscale grid cache: {len(grid_cache)} entries")
    print(f"  local fpocket cache: {len(local_cache)} entries")

    predictor = FpocketGridBoxPredictor()

    boxes: dict[str, dict] = {}
    n_cache_mid = n_cache_local = n_fpocket = n_centroid = n_skip = 0

    for i, pdb in enumerate(pdbs, 1):
        rid = pdb.stem
        if rid in grid_cache:
            e = grid_cache[rid]
            boxes[rid] = {"center": e["center"], "size": e["size"], "source": "fpocket-cache"}
            n_cache_mid += 1
            continue
        if rid in local_cache:
            boxes[rid] = local_cache[rid]
            n_cache_local += 1
            continue

        # fpocket 実行
        b = predict_fpocket(pdb, predictor)
        if b is not None:
            boxes[rid] = b
            local_cache[rid] = b
            n_fpocket += 1
        else:
            c = ca_centroid(pdb)
            if c is None:
                n_skip += 1
                continue
            boxes[rid] = {"center": list(c), "size": [CENTROID_SIZE] * 3, "source": "centroid"}
            n_centroid += 1

        if i % 50 == 0 or i == len(pdbs):
            print(f"  [{i:4d}/{len(pdbs)}] mid={n_cache_mid} local={n_cache_local} fp={n_fpocket} centroid={n_centroid} skip={n_skip}")
            CACHE_LOCAL.write_text(json.dumps(local_cache, indent=2))

    CACHE_LOCAL.write_text(json.dumps(local_cache, indent=2))
    BOXES_JSON.write_text(json.dumps(boxes, indent=2))
    with open(BOXES_TSV, "w") as f:
        f.write("id\tcx\tcy\tcz\tsx\tsy\tsz\n")
        for rid, b in boxes.items():
            f.write(f"{rid}\t" + "\t".join(map(str, list(b["center"]) + list(b["size"]))) + "\n")

    print(f"=== 完了 ===")
    print(f"  midscale cache:  {n_cache_mid}")
    print(f"  local fp cache:  {n_cache_local}")
    print(f"  fpocket new:     {n_fpocket}")
    print(f"  centroid (fall): {n_centroid}")
    print(f"  skipped:         {n_skip}")
    print(f"  boxes.tsv: {BOXES_TSV} ({len(boxes)} entries)")


if __name__ == "__main__":
    main()
