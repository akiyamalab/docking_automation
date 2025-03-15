#!/usr/bin/env python3
"""
Meekoを使用したリガンド準備とVinaドッキング計算のデモンストレーション

このスクリプトは、以下の処理を行います：
1. SDFファイルからリガンドを読み込む
2. Meekoを使ってリガンドをPDBQT形式に変換
3. AutoDock Vinaを使ってドッキング計算を実行
4. 結果を解析して表示
"""

import os
import sys
import time
import uuid
import logging
import tempfile
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple
import shutil

# RDKitとMeekoのインポート
from rdkit import Chem
from rdkit.Chem import AllChem
import rdkit.Chem.Descriptors  # 明示的にDescriptorsモジュールをインポート
from meeko import MoleculePreparation, PDBQTWriterLegacy
# Meekoの使用方法を修正

# ロギング設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("meeko_vina_docking")

# 現在の作業ディレクトリを取得
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

# 値オブジェクト
class MoleculeFormat:
    """分子形式を表す値オブジェクト"""
    def __init__(self, name: str):
        self.name = name

class GridBox:
    """ドッキング計算のグリッドボックスを表す値オブジェクト"""
    def __init__(self, center_x: float, center_y: float, center_z: float, 
                 size_x: float, size_y: float, size_z: float):
        self.center_x = center_x
        self.center_y = center_y
        self.center_z = center_z
        self.size_x = size_x
        self.size_y = size_y
        self.size_z = size_z
    
    def to_dict(self) -> dict:
        """辞書形式に変換"""
        return {
            'center_x': self.center_x,
            'center_y': self.center_y,
            'center_z': self.center_z,
            'size_x': self.size_x,
            'size_y': self.size_y,
            'size_z': self.size_z
        }

# エンティティ
class Molecule:
    """分子を表す基底エンティティ"""
    def __init__(self, id: str, path: str, format: MoleculeFormat):
        self.id = id
        self.path = path
        self.format = format
        self.properties = {}
    
    def get_path(self) -> str:
        """ファイルパスを取得"""
        return self.path

class Compound(Molecule):
    """化合物を表すエンティティ"""
    def __init__(self, id: str, path: str, format: MoleculeFormat, 
                 rdkit_mol: Optional[Chem.Mol] = None):
        super().__init__(id, path, format)
        self.rdkit_mol = rdkit_mol
        self.is_prepared = False
        self.molecular_weight = None
        
        # 分子量を計算（rdkit_molが提供されている場合）
        if rdkit_mol:
            self.calculate_properties()
    
    def calculate_properties(self):
        """RDKitを使用して分子特性を計算"""
        if self.rdkit_mol:
            self.molecular_weight = rdkit.Chem.Descriptors.MolWt(self.rdkit_mol)
            self.properties['molecular_weight'] = self.molecular_weight
            self.properties['logp'] = rdkit.Chem.Descriptors.MolLogP(self.rdkit_mol)
            self.properties['num_rotatable_bonds'] = rdkit.Chem.Descriptors.NumRotatableBonds(self.rdkit_mol)
            self.properties['num_hb_donors'] = rdkit.Chem.Descriptors.NumHDonors(self.rdkit_mol)
            self.properties['num_hb_acceptors'] = rdkit.Chem.Descriptors.NumHAcceptors(self.rdkit_mol)

class Protein(Molecule):
    """タンパク質を表すエンティティ"""
    def __init__(self, id: str, path: str, format: MoleculeFormat):
        super().__init__(id, path, format)
        self.is_prepared = format.name == "PDBQT"  # PDBQT形式なら準備済みとみなす

class Ligand:
    """ドッキング計算用のリガンドを表すエンティティ"""
    def __init__(self, compound: Compound, name: Optional[str] = None):
        self.compound = compound
        self.name = name or compound.id
    
    def is_prepared(self) -> bool:
        """リガンドが計算のために準備されているかを確認"""
        return self.compound.is_prepared and self.compound.format.name == "PDBQT"
    
    def get_path(self) -> str:
        """リガンドファイルのパス"""
        return self.compound.get_path()
    
    def get_molecular_weight(self) -> Optional[float]:
        """分子量を取得"""
        return self.compound.molecular_weight

class Receptor:
    """ドッキング計算用のレセプターを表すエンティティ"""
    def __init__(self, protein: Protein, name: Optional[str] = None):
        self.protein = protein
        self.name = name or protein.id
    
    def is_prepared(self) -> bool:
        """レセプターが計算のために準備されているかを確認"""
        return self.protein.is_prepared
    
    def get_path(self) -> str:
        """レセプターファイルのパス"""
        return self.protein.get_path()

class DockingTask:
    """ドッキング計算タスクを表すエンティティ"""
    def __init__(self, id: str, ligand: Ligand, receptor: Receptor, 
                 grid_box: GridBox, parameters: Dict[str, Any] = None):
        self.id = id
        self.ligand = ligand
        self.receptor = receptor
        self.grid_box = grid_box
        self.parameters = parameters or {}
        self.status = "PENDING"
        self.created_at = time.time()
        self.error_message = None
    
    def is_ready(self) -> bool:
        """タスクが実行可能かどうかを確認"""
        # 分子量フィルタリング（700以上はスキップ）
        mol_weight = self.ligand.get_molecular_weight()
        if mol_weight and mol_weight >= 700:
            self.error_message = f"Molecular weight too high: {mol_weight} >= 700"
            return False
            
        return (
            self.ligand.is_prepared() and 
            self.receptor.is_prepared()
        )
    
    def mark_as_running(self):
        """タスクを実行中にマーク"""
        self.status = "RUNNING"
    
    def mark_as_completed(self):
        """タスクを完了にマーク"""
        self.status = "COMPLETED"
    
    def mark_as_failed(self, error_message: str):
        """タスクを失敗にマーク"""
        self.status = "FAILED"
        self.error_message = error_message

class DockingResult:
    """ドッキング計算結果を表すエンティティ"""
    def __init__(self, id: str, task: DockingTask, output_path: str, log_path: str):
        self.id = id
        self.task = task
        self.output_path = output_path
        self.log_path = log_path
        self.poses = []
        self.scores = []
        self.execution_time = None
        self.created_at = time.time()
    
    def add_score(self, mode: int, affinity: float, rmsd_lb: float, rmsd_ub: float):
        """スコアを追加"""
        self.scores.append({
            'mode': mode,
            'affinity': affinity, 
            'rmsd_lb': rmsd_lb,
            'rmsd_ub': rmsd_ub
        })
    
    def get_best_score(self) -> Optional[Dict[str, Any]]:
        """最良のスコアを取得"""
        if not self.scores:
            return None
        return min(self.scores, key=lambda x: x['affinity'])

# サービス
class MeekoMoleculePreparationService:
    """Meekoを使用した分子準備サービス"""
    
    def prepare_ligand(self, compound: Compound) -> Compound:
        """化合物をリガンドとして準備
        
        Meekoを使用して、以下の処理を行います：
        1. 3D構造の最適化
        2. 回転可能な結合の検出
        3. 原子タイプの割り当て
        4. 部分電荷の計算
        5. PDBQTフォーマットへの変換
        """
        if compound.is_prepared:
            return compound
        
        # 出力ディレクトリの作成
        output_dir = os.path.join(CURRENT_DIR, "prepared_ligands")
        os.makedirs(output_dir, exist_ok=True)
        
        # 出力ファイルパスの設定
        output_path = os.path.join(output_dir, f"{compound.id}.pdbqt")
        
        try:
            # RDKitの分子オブジェクトが存在しない場合はファイルから読み込む
            if compound.rdkit_mol is None:
                if compound.format.name == "SDF":
                    mol = Chem.MolFromMolFile(compound.path)
                else:
                    logger.error(f"Unsupported format: {compound.format.name}")
                    raise ValueError(f"Unsupported format: {compound.format.name}")
                
                if mol is None:
                    logger.error(f"Failed to read molecule from {compound.path}")
                    raise ValueError(f"Failed to read molecule from {compound.path}")
                
                # 3D座標がない場合は生成
                if mol.GetNumConformers() == 0:
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol)
                    AllChem.UFFOptimizeMolecule(mol)
                
                compound.rdkit_mol = mol
                compound.calculate_properties()
            
            # Meekoによる準備
            # 明示的に水素原子を追加し、3D座標を持つことを確認
            mol = compound.rdkit_mol
            if mol.GetNumAtoms() > 0 and mol.GetNumConformers() > 0:
                # 明示的水素が必要
                if Chem.GetFormalCharge(mol) != sum([a.GetFormalCharge() for a in mol.GetAtoms()]):
                    mol = Chem.AddHs(mol)
                
                # Meekoの準備パラメータを設定
                preparator = MoleculePreparation()
                preparator.prepare(mol)
            else:
                raise ValueError("Molecule has no atoms or no 3D coordinates")
            
            # PDBQTファイルに書き出し
            # MoleculePreparationクラス自体がwrite_pdbqt_fileメソッドを持つ
            preparator.write_pdbqt_file(output_path)
            
            # 準備済みのCompoundオブジェクトを作成
            prepared_compound = Compound(
                id=compound.id,
                path=output_path,
                format=MoleculeFormat("PDBQT"),
                rdkit_mol=compound.rdkit_mol
            )
            prepared_compound.is_prepared = True
            prepared_compound.molecular_weight = compound.molecular_weight
            prepared_compound.properties = compound.properties.copy()
            
            return prepared_compound
            
        except Exception as e:
            logger.error(f"Error preparing ligand: {e}")
            raise
            
    def read_molecules_from_sdf(self, sdf_path: str, limit: int = None) -> List[Compound]:
        """SDFファイルから分子を読み込む"""
        logger.info(f"Reading molecules from {sdf_path}")
        
        compounds = []
        suppl = Chem.SDMolSupplier(sdf_path)
        
        # limit数まで読み込み、無効な分子はスキップ
        count = 0
        for mol_idx, mol in enumerate(suppl):
            if mol is not None:
                # 明示的な水素原子を追加し、3D座標を生成
                try:
                    # 明示的に水素原子を追加
                    mol = Chem.AddHs(mol)
                    
                    # 3D座標を生成
                    embed_result = AllChem.EmbedMolecule(mol, randomSeed=42)
                    if embed_result == -1:
                        logger.warning(f"Failed to generate 3D coordinates for molecule {mol_idx} with standard parameters")
                        # 代替パラメータで再試行
                        params = AllChem.ETKDGv2()
                        params.useRandomCoords = True
                        AllChem.EmbedMolecule(mol, params)
                    
                    # 分子構造を最適化
                    AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                except Exception as e:
                    logger.warning(f"Failed to generate 3D coordinates for molecule {mol_idx}: {e}")
                    continue
                
                # 分子名または化合物IDを取得
                name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Compound_{mol_idx}"
                
                # 一時ファイルに書き出し
                output_dir = os.path.join(CURRENT_DIR, "input_ligands")
                os.makedirs(output_dir, exist_ok=True)
                mol_path = os.path.join(output_dir, f"{name}.sdf")
                
                # SDFに保存
                writer = Chem.SDWriter(mol_path)
                writer.write(mol)
                writer.close()
                
                # Compoundオブジェクトを作成
                compound = Compound(
                    id=name,
                    path=mol_path,
                    format=MoleculeFormat("SDF"),
                    rdkit_mol=mol
                )
                
                compounds.append(compound)
                count += 1
                
                if limit and count >= limit:
                    break
        
        logger.info(f"Read {len(compounds)} valid molecules")
        return compounds

class VinaDockingService:
    """AutoDock Vinaを使用したドッキング計算サービス"""
    
    def __init__(self, vina_path: str = "vina"):
        self.vina_path = vina_path
    
    def execute(self, task: DockingTask) -> DockingResult:
        """ドッキング計算を実行"""
        # タスクの妥当性を検証
        if not task.is_ready():
            error_msg = task.error_message or "Task is not ready for execution"
            logger.error(error_msg)
            task.mark_as_failed(error_msg)
            raise ValueError(error_msg)
        
        # 出力ディレクトリの作成
        output_dir = os.path.join(CURRENT_DIR, "docking_results")
        os.makedirs(output_dir, exist_ok=True)
        
        # 出力ファイルのパス
        output_path = os.path.join(output_dir, f"result_{task.id}.pdbqt")
        log_path = os.path.join(output_dir, f"log_{task.id}.txt")
        
        # タスクを実行中にマーク
        task.mark_as_running()
        
        try:
            start_time = time.time()
            
            # 一時設定ファイルの作成
            config_path = self._create_config_file(task.grid_box, task.parameters)
            
            # Vinaコマンドの構築
            cmd = [
                self.vina_path,
                "--receptor", task.receptor.get_path(),
                "--ligand", task.ligand.get_path(),
                "--config", config_path,
                "--out", output_path,
                "--log", log_path
            ]
            
            # CPUパラメータの追加（存在する場合）
            if 'cpu' in task.parameters:
                cmd.extend(["--cpu", str(task.parameters['cpu'])])
            
            # コマンドを実行
            logger.info(f"Running Vina command: {' '.join(cmd)}")
            
            process = subprocess.run(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            
            # 処理時間を計算
            end_time = time.time()
            execution_time = end_time - start_time
            
            # 結果を確認
            if process.returncode != 0:
                error_msg = f"Vina execution failed: {process.stderr}"
                logger.error(error_msg)
                task.mark_as_failed(error_msg)
                raise RuntimeError(error_msg)
            
            logger.info(f"Vina calculation completed in {execution_time:.2f} seconds")
            
            # 出力ファイルの存在確認
            if not os.path.exists(output_path):
                error_msg = "Output file was not created"
                logger.error(error_msg)
                task.mark_as_failed(error_msg)
                raise RuntimeError(error_msg)
            
            # タスクを完了にマーク
            task.mark_as_completed()
            
            # 結果オブジェクトの作成
            result = DockingResult(
                id=str(uuid.uuid4()),
                task=task,
                output_path=output_path,
                log_path=log_path
            )
            result.execution_time = execution_time
            
            # ログファイルからスコアを解析
            self._parse_scores(result)
            
            return result
            
        except Exception as e:
            logger.error(f"Error executing docking task: {e}")
            task.mark_as_failed(str(e))
            raise
            
        finally:
            # 一時ファイルの削除
            if 'config_path' in locals() and os.path.exists(config_path):
                os.unlink(config_path)
    
    def create_task(self, ligand: Ligand, receptor: Receptor, 
                    grid_box: GridBox, parameters: Dict[str, Any] = None) -> DockingTask:
        """ドッキングタスクを作成"""
        # リガンドとレセプターの妥当性を確認
        if not ligand.is_prepared():
            raise ValueError("Ligand is not properly prepared for docking")
        
        if not receptor.is_prepared():
            raise ValueError("Receptor is not properly prepared for docking")
        
        # タスクの作成
        task = DockingTask(
            id=str(uuid.uuid4())[:8],
            ligand=ligand,
            receptor=receptor,
            grid_box=grid_box,
            parameters=parameters or {}
        )
        
        return task
    
    def _create_config_file(self, grid_box: GridBox, parameters: Dict[str, Any]) -> str:
        """Vinaの設定ファイルを作成"""
        # 一時ファイルを作成
        fd, config_path = tempfile.mkstemp(suffix=".conf", prefix="vina_")
        with os.fdopen(fd, 'w') as f:
            # グリッドボックスの設定
            f.write(f"center_x = {grid_box.center_x}\n")
            f.write(f"center_y = {grid_box.center_y}\n")
            f.write(f"center_z = {grid_box.center_z}\n")
            f.write(f"size_x = {grid_box.size_x}\n")
            f.write(f"size_y = {grid_box.size_y}\n")
            f.write(f"size_z = {grid_box.size_z}\n")
            
            # その他のパラメータ
            vina_params = {
                'exhaustiveness': parameters.get('exhaustiveness', 8),
                'num_modes': parameters.get('num_modes', 9),
                'energy_range': parameters.get('energy_range', 3)
            }
            
            # コマンドライン引数として渡すパラメータは設定ファイルに書き込まない
            skip_params = ["cpu"]
            
            for name, value in vina_params.items():
                if name not in skip_params:
                    f.write(f"{name} = {value}\n")
        
        return config_path
    
    def _parse_scores(self, result: DockingResult) -> None:
        """ログファイルからスコアを解析"""
        try:
            with open(result.log_path, 'r') as f:
                log_content = f.read()
                
                # スコアテーブルの位置を特定
                import re
                table_start = log_content.find("-----+------------+----------+----------")
                if table_start == -1:
                    logger.warning("Score table not found in log file")
                    return
                
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
                        
                        result.add_score(mode_num, affinity, rmsd_lb, rmsd_ub)
        
        except Exception as e:
            logger.error(f"Error parsing log file: {e}")

def main():
    """メイン関数"""
    # 入力ファイルの確認
    sdf_path = os.path.join(CURRENT_DIR, "sample_compounds.sdf")
    receptor_path = os.path.join(CURRENT_DIR, "protein.pdbqt")
    
    if not os.path.exists(sdf_path):
        logger.error(f"SDF file not found: {sdf_path}")
        return 1
    
    if not os.path.exists(receptor_path):
        logger.error(f"Receptor file not found: {receptor_path}")
        return 1
    
    # 結果ディレクトリの作成
    output_dir = os.path.join(CURRENT_DIR, "meeko_vina_results")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 1. サービスの初期化
        preparation_service = MeekoMoleculePreparationService()
        docking_service = VinaDockingService(vina_path="/advina/1.1.2/bin/vina")  # Vinaの実行パス
        
        # 2. レセプターの準備
        logger.info(f"Loading receptor: {receptor_path}")
        receptor_protein = Protein(
            id="receptor",
            path=receptor_path,
            format=MoleculeFormat("PDBQT")
        )
        receptor = Receptor(protein=receptor_protein)
        
        # 3. SDFファイルから化合物を読み込む
        compounds = preparation_service.read_molecules_from_sdf(sdf_path, limit=3)  # 最初の3つだけ処理
        
        # 4. ドッキング設定の準備
        grid_box = GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
        
        parameters = {
            'exhaustiveness': 8,
            'num_modes': 5,
            'energy_range': 3,
            'cpu': 1
        }
        
        # 5. 各化合物に対してドッキング計算を実行
        results = []
        skipped = 0
        
        for compound in compounds:
            logger.info(f"Processing compound: {compound.id}")
            
            # 分子量の確認（700以上はスキップ）
            if compound.molecular_weight and compound.molecular_weight >= 700:
                logger.info(f"Skipping compound {compound.id} due to high molecular weight: {compound.molecular_weight:.2f}")
                skipped += 1
                continue
            
            # リガンド準備
            try:
                logger.info("Preparing ligand with Meeko")
                prepared_compound = preparation_service.prepare_ligand(compound)
                ligand = Ligand(compound=prepared_compound)
                
                # ドッキングタスクの作成
                task = docking_service.create_task(
                    ligand=ligand,
                    receptor=receptor,
                    grid_box=grid_box,
                    parameters=parameters
                )
                
                # ドッキング計算の実行
                logger.info(f"Executing docking task: {task.id}")
                result = docking_service.execute(task)
                
                # 結果を保存
                results.append(result)
                
                # 結果の表示
                logger.info(f"Docking completed for {compound.id}")
                best_score = result.get_best_score()
                if best_score:
                    logger.info(f"Best score: {best_score['affinity']:.2f} kcal/mol")
                logger.info(f"Execution time: {result.execution_time:.2f} seconds")
                
            except Exception as e:
                logger.error(f"Error processing compound {compound.id}: {e}")
                continue
        
        # 6. 結果の要約を保存
        summary_path = os.path.join(output_dir, "docking_summary.txt")
        with open(summary_path, "w") as f:
            f.write(f"Docking Results Summary\n")
            f.write(f"======================\n\n")
            f.write(f"Total compounds processed: {len(compounds)}\n")
            f.write(f"Successful dockings: {len(results)}\n")
            f.write(f"Skipped compounds (molecular weight >= 700): {skipped}\n\n")
            
            f.write(f"Best scores:\n")
            f.write(f"-----------\n")
            
            for result in sorted(results, key=lambda r: r.get_best_score()['affinity'] if r.get_best_score() else 0):
                best_score = result.get_best_score()
                if best_score:
                    compound_name = result.task.ligand.name
                    f.write(f"{compound_name}: {best_score['affinity']:.2f} kcal/mol\n")
            
            f.write(f"\nDetailed results:\n")
            f.write(f"----------------\n")
            
            for result in results:
                compound_name = result.task.ligand.name
                f.write(f"\nCompound: {compound_name}\n")
                f.write(f"Molecular Weight: {result.task.ligand.get_molecular_weight():.2f}\n")
                f.write(f"Execution Time: {result.execution_time:.2f} seconds\n")
                f.write(f"Output File: {os.path.basename(result.output_path)}\n")
                f.write(f"Scores:\n")
                
                for score in result.scores:
                    f.write(f"  Mode {score['mode']}: {score['affinity']:.2f} kcal/mol (RMSD: {score['rmsd_lb']:.2f})\n")
        
        logger.info(f"Results summary saved to: {summary_path}")
        return 0
        
    except Exception as e:
        logger.error(f"Error in main process: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())