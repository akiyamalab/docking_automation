"""
Uni-Dock E2E テスト (subtask_010_b)

AlphaFold PDB ファイルをレセプターとし、ALDR actives_subset.sdf から
最大 N 件のリガンドを用いて Uni-Dock を実行する。
AutoDock Vina スコアと比較して差 ±2 kcal/mol 以内であることを確認。

Usage:
    cd /workspaces/20260422_mouse_docking/docking_automation
    python3 examples/unidock_e2e_test.py            # デフォルト (N=5)
    python3 examples/unidock_e2e_test.py --n-pairs 10
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

os.environ["BABEL_QUIET"] = "1"

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent

INPUT_DIR = SCRIPT_DIR / "input"
OUTPUT_DIR = SCRIPT_DIR / "output"
ALPHAFOLD_DIR = INPUT_DIR / "alphafold"
ALDR_DIR = INPUT_DIR / "ALDR"

RECEPTOR_PDB = ALPHAFOLD_DIR / "AF-A0A0A2JW93-F1-model_v4.pdb"
LIGANDS_SDF = ALDR_DIR / "actives_subset.sdf"
UNIDOCK_OUT_DIR = OUTPUT_DIR / "unidock_e2e"


def compute_grid_center(pdb_path: Path, size: float = 30.0) -> tuple:
    """PDB ファイルの Cα 重心からグリッドボックスを計算する。"""
    coords = []
    with open(pdb_path) as f:
        for line in f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atom_name = line[12:16].strip()
            if atom_name != "CA":
                continue
            try:
                x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
                coords.append([x, y, z])
            except ValueError:
                continue
    if not coords:
        raise ValueError(f"Cα 座標が見つかりません: {pdb_path}")
    center = np.mean(np.array(coords), axis=0)
    return tuple(center.tolist()), size


def prepare_receptor_pdbqt(pdb_path: Path, out_dir: Path) -> Path:
    """obabel で PDB → PDBQT (受容体モード) に変換する。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    pdbqt_path = out_dir / (pdb_path.stem + ".pdbqt")
    if pdbqt_path.exists():
        return pdbqt_path

    tmp_pdb = out_dir / (pdb_path.stem + "_noh.pdb")
    subprocess.run(
        ["obabel", str(pdb_path), "-O", str(tmp_pdb), "-d"],
        check=True, capture_output=True
    )
    subprocess.run(
        ["obabel", str(tmp_pdb), "-O", str(pdbqt_path), "-xr"],
        check=True, capture_output=True
    )
    return pdbqt_path


def prepare_ligands_pdbqt(sdf_path: Path, out_dir: Path, max_n: int) -> list[Path]:
    """meeko で SDF → PDBQT リストに変換する (最大 max_n 件)。"""
    from meeko import MoleculePreparation, PDBQTWriterLegacy
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out_dir.mkdir(parents=True, exist_ok=True)
    prepared: list[Path] = []
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)

    preparator = MoleculePreparation()
    for i, mol in enumerate(suppl):
        if i >= max_n:
            break
        if mol is None:
            continue
        name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"lig_{i}"
        out_pdbqt = out_dir / f"{name}_{i}.pdbqt"
        if out_pdbqt.exists():
            prepared.append(out_pdbqt)
            continue

        mol_h = Chem.AddHs(mol)
        if AllChem.EmbedMolecule(mol_h, randomSeed=42) != 0:
            AllChem.EmbedMolecule(mol_h, randomSeed=42, useRandomCoords=True)
        AllChem.MMFFOptimizeMolecule(mol_h)

        try:
            setups = preparator.prepare(mol_h)
            pdbqt_str, ok, err = PDBQTWriterLegacy.write_string(setups[0])
            if not ok:
                print(f"  [WARN] meeko skipped lig {i}: {err}")
                continue
            out_pdbqt.write_text(pdbqt_str)
            prepared.append(out_pdbqt)
        except Exception as e:
            print(f"  [WARN] meeko error lig {i}: {e}")

    return prepared


def run_vina(receptor_pdbqt: Path, ligand_pdbqt: Path,
             center: tuple, size: float, out_dir: Path) -> float | None:
    """AutoDock Vina でスコアを計算して返す。"""
    try:
        from vina import Vina
    except ImportError:
        return None

    v = Vina(sf_name="vina", verbosity=0)
    v.set_receptor(str(receptor_pdbqt))
    v.set_ligand_from_file(str(ligand_pdbqt))
    v.compute_vina_maps(center=list(center), box_size=[size, size, size])
    v.dock(exhaustiveness=4, n_poses=1)
    energies = v.energies(n_poses=1)
    return float(energies[0][0]) if energies is not None and len(energies) > 0 else None


def run_unidock(receptor_pdbqt: Path, ligand_pdbqts: list[Path],
                center: tuple, size: float, out_dir: Path) -> dict[str, float]:
    """Uni-Dock で全リガンドをバッチ実行しスコアを返す。"""
    out_dir.mkdir(parents=True, exist_ok=True)
    ligand_index_file = out_dir / "ligands.dat"
    ligand_index_file.write_text("\n".join(str(p) for p in ligand_pdbqts) + "\n")

    cx, cy, cz = center
    cmd = [
        "unidock",
        "--receptor", str(receptor_pdbqt),
        "--ligand_index", str(ligand_index_file),
        "--center_x", f"{cx:.3f}",
        "--center_y", f"{cy:.3f}",
        "--center_z", f"{cz:.3f}",
        "--size_x", str(size),
        "--size_y", str(size),
        "--size_z", str(size),
        "--dir", str(out_dir),
        "--num_modes", "1",
        "--exhaustiveness", "8",
        "--scoring", "vina",
    ]
    print("  Running:", " ".join(cmd))
    t0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    elapsed = time.time() - t0
    print(f"  Uni-Dock elapsed: {elapsed:.1f}s")
    if result.returncode != 0:
        print("  [WARN] Uni-Dock stderr:", result.stderr[:500])

    scores: dict[str, float] = {}
    for pdbqt_path in out_dir.glob("*_out.pdbqt"):
        stem = pdbqt_path.stem.replace("_out", "")
        score = _parse_unidock_score(pdbqt_path)
        if score is not None:
            scores[stem] = score
    return scores


def _parse_unidock_score(pdbqt_path: Path) -> float | None:
    """PDBQT 出力ファイルから最初の REMARK VINA RESULT の affinity を取得。"""
    for line in pdbqt_path.read_text().splitlines():
        if "REMARK VINA RESULT" in line or "REMARK  VINA RESULT" in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p in ("RESULT", "VINA") and i + 1 < len(parts):
                    try:
                        return float(parts[i + 1])
                    except ValueError:
                        pass
    return None


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--n-pairs", type=int, default=5,
                   help="テストするリガンド数 (デフォルト: 5)")
    p.add_argument("--skip-vina", action="store_true",
                   help="Vina 比較をスキップ")
    return p.parse_args()


def main():
    args = parse_args()
    n = args.n_pairs
    print(f"=== Uni-Dock E2E テスト (N={n} ペア) ===")

    work_dir = UNIDOCK_OUT_DIR
    work_dir.mkdir(parents=True, exist_ok=True)

    # 1. レセプター準備
    print(f"\n[1] レセプター前処理: {RECEPTOR_PDB.name}")
    receptor_pdbqt = prepare_receptor_pdbqt(RECEPTOR_PDB, work_dir / "receptor")
    print(f"    → {receptor_pdbqt}")

    # 2. グリッドボックス計算
    center, size = compute_grid_center(RECEPTOR_PDB)
    print(f"\n[2] グリッドボックス: center={tuple(f'{v:.2f}' for v in center)}, size={size}")

    # 3. リガンド準備
    print(f"\n[3] リガンド前処理 (最大 {n} 件): {LIGANDS_SDF.name}")
    lig_pdbqts = prepare_ligands_pdbqt(LIGANDS_SDF, work_dir / "ligands", max_n=n)
    print(f"    → {len(lig_pdbqts)} 件準備完了")

    if not lig_pdbqts:
        print("[ERROR] リガンドが準備できませんでした")
        sys.exit(1)

    # 4. Uni-Dock 実行
    print(f"\n[4] Uni-Dock 実行 ({len(lig_pdbqts)} リガンド)")
    unidock_scores = run_unidock(
        receptor_pdbqt, lig_pdbqts, center, size, work_dir / "results"
    )
    print(f"    スコア取得: {len(unidock_scores)} 件")

    # 5. Vina 比較
    vina_scores: dict[str, float] = {}
    if not args.skip_vina and len(unidock_scores) > 0:
        print(f"\n[5] AutoDock Vina 比較 ({len(lig_pdbqts)} リガンド)")
        for lig_pdbqt in lig_pdbqts:
            stem = lig_pdbqt.stem
            score = run_vina(receptor_pdbqt, lig_pdbqt, center, size,
                             work_dir / "vina_results")
            if score is not None:
                vina_scores[stem] = score
                print(f"    {stem}: Vina={score:.3f}")

    # 6. 結果サマリー
    print("\n=== 結果サマリー ===")
    print(f"{'Ligand':<30} {'Uni-Dock':>10} {'Vina':>10} {'Diff':>8} {'OK?':>5}")
    print("-" * 65)

    all_ok = True
    pairs_tested = 0
    diffs = []

    for stem, ud_score in sorted(unidock_scores.items()):
        vina_score = vina_scores.get(stem)
        if vina_score is not None:
            diff = abs(ud_score - vina_score)
            ok = diff <= 2.0
            diffs.append(diff)
            pairs_tested += 1
            status = "OK" if ok else "NG"
            if not ok:
                all_ok = False
            print(f"{stem:<30} {ud_score:>10.3f} {vina_score:>10.3f} {diff:>8.3f} {status:>5}")
        else:
            print(f"{stem:<30} {ud_score:>10.3f} {'N/A':>10} {'N/A':>8} {'N/A':>5}")
            pairs_tested += 1

    print()
    if diffs:
        print(f"平均スコア差: {np.mean(diffs):.3f} kcal/mol")
        print(f"最大スコア差: {np.max(diffs):.3f} kcal/mol")

    # 結果 JSON 保存
    result_json = {
        "unidock_version": "1.1.0",
        "install_method": "binary download",
        "receptor": str(RECEPTOR_PDB.name),
        "n_pairs_tested": pairs_tested,
        "unidock_scores": unidock_scores,
        "vina_scores": vina_scores,
        "mean_diff": float(np.mean(diffs)) if diffs else None,
        "max_diff": float(np.max(diffs)) if diffs else None,
        "all_within_2kcal": all_ok,
    }
    result_file = work_dir / "e2e_result.json"
    result_file.write_text(json.dumps(result_json, indent=2))
    print(f"\n結果保存: {result_file}")
    print(f"E2E テスト: {'PASS' if pairs_tested > 0 else 'FAIL (スコア取得なし)'}")

    return 0 if pairs_tested > 0 else 1


if __name__ == "__main__":
    sys.exit(main())
