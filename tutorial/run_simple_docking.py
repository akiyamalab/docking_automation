#!/usr/bin/env python3
"""
単純な分子を使ったドッキング計算のデモンストレーション

1. RDKitを使って単純な分子を生成
2. Meekoでリガンド準備
3. AutoDock Vinaでドッキング計算
"""

import os
import sys
import time
import uuid
import logging
import tempfile
import subprocess
from pathlib import Path
import shutil

# RDKitとMeekoのインポート
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from meeko import MoleculePreparation, PDBQTWriterLegacy

# ロギング設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("simple_docking")

# 現在の作業ディレクトリを取得
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def create_simple_molecule(smiles: str, mol_id: str):
    """単純な分子を作成してPDB形式に変換"""
    output_dir = os.path.join(CURRENT_DIR, "input_ligands")
    os.makedirs(output_dir, exist_ok=True)
    
    # SMILESから分子を作成
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logger.error(f"Failed to create molecule from SMILES: {smiles}")
        return None
    
    # 分子名を設定
    mol.SetProp("_Name", mol_id)
    
    # 明示的な水素原子を追加
    mol = Chem.AddHs(mol)
    
    # 3D座標を生成
    try:
        embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
        if embed_result == -1:
            logger.warning("Failed to generate 3D coordinates with standard parameters")
            # 代替パラメータで再試行
            params = AllChem.ETKDGv2()
            params.useRandomCoords = True
            AllChem.EmbedMolecule(mol, params)
        
        # 分子構造を最適化
        AllChem.UFFOptimizeMolecule(mol, maxIters=500)
    except Exception as e:
        logger.error(f"Error generating 3D coordinates: {e}")
        return None
    
    # 分子特性を出力
    logger.info(f"Molecular weight: {Descriptors.MolWt(mol):.2f}")
    logger.info(f"LogP: {Descriptors.MolLogP(mol):.2f}")
    logger.info(f"Rotatable bonds: {Descriptors.NumRotatableBonds(mol)}")
    logger.info(f"H-bond donors: {Descriptors.NumHDonors(mol)}")
    logger.info(f"H-bond acceptors: {Descriptors.NumHAcceptors(mol)}")
    
    # PDBファイルに出力
    pdb_path = os.path.join(output_dir, f"{mol_id}.pdb")
    Chem.MolToPDBFile(mol, pdb_path)
    logger.info(f"Molecule saved to: {pdb_path}")
    
    return {
        "mol": mol,
        "path": pdb_path,
        "id": mol_id,
        "properties": {
            "molecular_weight": Descriptors.MolWt(mol),
            "logp": Descriptors.MolLogP(mol),
            "num_rotatable_bonds": Descriptors.NumRotatableBonds(mol),
            "num_hb_donors": Descriptors.NumHDonors(mol),
            "num_hb_acceptors": Descriptors.NumHAcceptors(mol)
        }
    }

def prepare_ligand_with_meeko(molecule_data):
    """Meekoを使ってリガンドを準備"""
    output_dir = os.path.join(CURRENT_DIR, "prepared_ligands")
    os.makedirs(output_dir, exist_ok=True)
    
    mol = molecule_data["mol"]
    mol_id = molecule_data["id"]
    output_path = os.path.join(output_dir, f"{mol_id}.pdbqt")
    
    try:
        # Meekoの準備
        preparator = MoleculePreparation()
        preparator.prepare(mol)
        
        # PDBQTファイルに書き出し - Meekoの書き出し方法
        preparator.write_pdbqt_file(output_path)
        
        logger.info(f"Prepared ligand saved to: {output_path}")
        
        return {
            "path": output_path,
            "id": mol_id,
            "is_prepared": True,
            "properties": molecule_data["properties"]
        }
    except Exception as e:
        logger.error(f"Error preparing ligand: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None

def run_vina_docking(ligand_path, receptor_path, output_dir, center_x=0, center_y=0, center_z=0, 
                     size_x=20, size_y=20, size_z=20, exhaustiveness=8):
    """AutoDock Vinaを使ってドッキング計算を実行"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 出力ファイルのパス
    output_path = os.path.join(output_dir, "docking_result.pdbqt")
    log_path = os.path.join(output_dir, "docking_log.txt")
    
    # 一時設定ファイルの作成
    fd, config_path = tempfile.mkstemp(suffix=".conf", prefix="vina_")
    with os.fdopen(fd, 'w') as f:
        # グリッドボックスの設定
        f.write(f"center_x = {center_x}\n")
        f.write(f"center_y = {center_y}\n")
        f.write(f"center_z = {center_z}\n")
        f.write(f"size_x = {size_x}\n")
        f.write(f"size_y = {size_y}\n")
        f.write(f"size_z = {size_z}\n")
        
        # その他のパラメータ
        f.write(f"exhaustiveness = {exhaustiveness}\n")
        f.write("num_modes = 9\n")
        f.write("energy_range = 3\n")
    
    # Vinaパスの設定
    vina_path = "/advina/1.1.2/bin/vina"
    
    # Vinaコマンドの構築
    cmd = [
        vina_path,
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--config", config_path,
        "--out", output_path,
        "--log", log_path,
        "--cpu", "1"
    ]
    
    try:
        # コマンドを実行
        logger.info(f"Running Vina command: {' '.join(cmd)}")
        start_time = time.time()
        
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # 結果を確認
        if process.returncode != 0:
            logger.error(f"Vina execution failed: {process.stderr}")
            return None
        
        logger.info(f"Vina calculation completed in {execution_time:.2f} seconds")
        
        # 結果を解析
        scores = parse_vina_log(log_path)
        
        return {
            "output_path": output_path,
            "log_path": log_path,
            "execution_time": execution_time,
            "scores": scores
        }
    
    except Exception as e:
        logger.error(f"Error executing Vina: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None
    
    finally:
        # 一時ファイルの削除
        if os.path.exists(config_path):
            os.unlink(config_path)

def parse_vina_log(log_path):
    """Vinaのログファイルを解析してスコアを抽出"""
    scores = []
    
    try:
        with open(log_path, 'r') as f:
            log_content = f.read()
            
            # スコアテーブルの位置を特定
            import re
            table_start = log_content.find("-----+------------+----------+----------")
            if table_start == -1:
                logger.warning("Score table not found in log file")
                return scores
            
            # モードのスコアを抽出
            # 例: "   1       -8.7      0.000      0.000"
            mode_pattern = r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)"
            lines = log_content[table_start:].split('\n')
            
            for line in lines[1:]:  # ヘッダー行をスキップ
                match = re.match(mode_pattern, line.strip())
                if match:
                    mode_num = int(match.group(1))
                    affinity = float(match.group(2))
                    rmsd_lb = float(match.group(3))
                    rmsd_ub = float(match.group(4))
                    
                    scores.append({
                        'mode': mode_num,
                        'affinity': affinity,
                        'rmsd_lb': rmsd_lb,
                        'rmsd_ub': rmsd_ub
                    })
    
    except Exception as e:
        logger.error(f"Error parsing log file: {e}")
    
    return scores

def main():
    # 出力ディレクトリの作成
    output_dir = os.path.join(CURRENT_DIR, "simple_docking_results")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # レセプターファイルの確認
    receptor_path = os.path.join(CURRENT_DIR, "protein.pdbqt")
    if not os.path.exists(receptor_path):
        logger.error(f"Receptor file not found: {receptor_path}")
        return 1
    
    # サンプル分子の作成
    # アスピリン（薬の代表例）
    aspirin_smiles = "CC(=O)OC1=CC=CC=C1C(=O)O"
    aspirin_data = create_simple_molecule(aspirin_smiles, "aspirin")
    
    # イブプロフェン（抗炎症薬）
    ibuprofen_smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    ibuprofen_data = create_simple_molecule(ibuprofen_smiles, "ibuprofen")
    
    # パラセタモール（鎮痛薬）
    paracetamol_smiles = "CC(=O)NC1=CC=C(C=C1)O"
    paracetamol_data = create_simple_molecule(paracetamol_smiles, "paracetamol")
    
    molecules = [m for m in [aspirin_data, ibuprofen_data, paracetamol_data] if m is not None]
    if not molecules:
        logger.error("No valid molecules could be created")
        return 1
    
    # 結果のまとめ
    results = []
    
    # 各分子でドッキング計算を実行
    for molecule_data in molecules:
        logger.info(f"Processing {molecule_data['id']}...")
        
        # リガンド準備
        prepared_ligand = prepare_ligand_with_meeko(molecule_data)
        if not prepared_ligand:
            logger.error(f"Failed to prepare {molecule_data['id']}")
            continue
        
        # ドッキング計算
        molecule_output_dir = os.path.join(output_dir, molecule_data['id'])
        docking_result = run_vina_docking(
            ligand_path=prepared_ligand["path"],
            receptor_path=receptor_path,
            output_dir=molecule_output_dir
        )
        
        if not docking_result:
            logger.error(f"Docking failed for {molecule_data['id']}")
            continue
        
        # 結果を保存
        results.append({
            "id": molecule_data['id'],
            "properties": molecule_data['properties'],
            "docking_result": docking_result
        })
    
    # 結果の要約を保存
    summary_path = os.path.join(output_dir, "docking_summary.txt")
    with open(summary_path, "w") as f:
        f.write(f"Docking Results Summary\n")
        f.write(f"======================\n\n")
        f.write(f"Total compounds processed: {len(molecules)}\n")
        f.write(f"Successful dockings: {len(results)}\n\n")
        
        f.write(f"Best scores:\n")
        f.write(f"-----------\n")
        
        for result in sorted(results, key=lambda r: r["docking_result"]["scores"][0]["affinity"] if r["docking_result"]["scores"] else 0):
            if result["docking_result"]["scores"]:
                best_score = result["docking_result"]["scores"][0]
                f.write(f"{result['id']}: {best_score['affinity']:.2f} kcal/mol\n")
        
        f.write(f"\nDetailed results:\n")
        f.write(f"----------------\n")
        
        for result in results:
            f.write(f"\nCompound: {result['id']}\n")
            f.write(f"Molecular Weight: {result['properties']['molecular_weight']:.2f}\n")
            f.write(f"LogP: {result['properties']['logp']:.2f}\n")
            f.write(f"Rotatable Bonds: {result['properties']['num_rotatable_bonds']}\n")
            f.write(f"Execution Time: {result['docking_result']['execution_time']:.2f} seconds\n")
            f.write(f"Scores:\n")
            
            for score in result["docking_result"]["scores"]:
                f.write(f"  Mode {score['mode']}: {score['affinity']:.2f} kcal/mol (RMSD: {score['rmsd_lb']:.2f})\n")
    
    logger.info(f"Results summary saved to: {summary_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())