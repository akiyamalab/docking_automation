#!/usr/bin/env python
# coding: utf-8

"""
ALDRタンパク質のドッキング計算を実行するスクリプト（自動グリッドボックス計算版）

分子構造データからドッキング中心とボックスサイズを自動計算して使用します。
2つの方法を示します：
1. 結晶リガンド（MOL2ファイル）から計算する方法
2. レセプター（タンパク質）構造から計算する方法
"""

import os
import sys
from pathlib import Path
import random
import uuid
import time
from rdkit import Chem
from rdkit.Chem import AllChem

from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.entity.compound import Compound
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.docking_task import DockingTask
from docking_automation.docking.entity.docking_result import DockingResult
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.docking_parameter import DockingParameter, ParameterType
from docking_automation.docking.value_object.docking_parameters import DockingParameters
from docking_automation.docking.value_object.pose import Pose
from docking_automation.docking.value_object.score import Score, ScoreType
from docking_automation.docking.service.docking_service import DockingService
from docking_automation.molecule.service.grid_box_calculator import GridBoxCalculator
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond


# 型変換は行わず、直接パラメータを渡す方式に変更
# 本番コードではなくチュートリアル実装として、型チェックエラーは無視します

# MOL2ファイルから原子座標を読み取り、MoleculeStructureを作成する関数
def load_mol2_structure(filename):
    """MOL2ファイルから分子構造を読み込む
    
    この簡易実装では、MOL2ファイルの原子座標のみを読み取ります。
    実際の実装ではRDKitやmeekoを使用することを推奨します。
    
    Args:
        filename: MOL2ファイルのパス
        
    Returns:
        MoleculeStructure: 分子構造オブジェクト
    """
    atoms = []
    atom_section = False
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                
                if line == '@<TRIPOS>ATOM':
                    atom_section = True
                    continue
                elif line.startswith('@<TRIPOS>'):
                    atom_section = False
                    continue
                
                if atom_section and line:
                    parts = line.split()
                    if len(parts) >= 5:  # 少なくとも原子ID、x、y、z座標が必要
                        atom_id = int(parts[0])
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                        element = parts[5].split('.')[0] if '.' in parts[5] else parts[5]
                        
                        atoms.append(Atom(
                            atom_id=atom_id,
                            element=element,
                            x=x, y=y, z=z,
                            charge=0.0,
                            properties={}
                        ))
        
        if not atoms:
            raise ValueError(f"エラー: {filename}から原子を読み取れませんでした")
        
        # 結合情報なしでMoleculeStructureを作成
        return MoleculeStructure(atoms=atoms, bonds=[])
    
    except Exception as e:
        print(f"MOL2ファイルの読み込み中にエラーが発生しました: {e}")
        raise


# モック実装を作成（本来はVinaのサービスを使用するが、ここではモックで代用）
class MockDockingService(DockingService):
    """デモ用のモックドッキングサービス"""
    
    def execute(self, task: DockingTask) -> DockingResult:
        """ドッキング計算を実行（モック）"""
        print(f"[モック] {task.ligand.name}のドッキング計算を実行中...")
        # ランダムなスコアを生成（既存コードに合わせたパラメータで）
        scores = []
        for i in range(3):
            # Scoreクラスの正しいパラメータを使用
            score = Score(
                value=random.uniform(-12.0, -5.0),
                score_type=ScoreType.BINDING_AFFINITY,  # 正しい列挙型を使用
                name=f"Pose_{i+1}"  # 名前は必須
            )
            scores.append(score)
        
        # モックのポーズを生成（既存コードに合わせたパラメータで）
        center = task.configuration.grid_box.get_center()
        poses = []
        for i in range(3):
            # モックの構造を作成（中心座標に配置）
            atom = Atom(
                atom_id=1,
                element="C",
                x=center[0],
                y=center[1],
                z=center[2]
            )
            # 構造を作成
            mock_structure = MoleculeStructure(atoms=[atom], bonds=[])
            
            # ポーズを作成（rankパラメータは必須）
            poses.append(Pose(structure=mock_structure, rank=i+1))
        
        # 結果を返す（既存コードに合わせたパラメータで）
        result = DockingResult(
            id=str(uuid.uuid4()),
            task=task,
            created_at=time.time(),
            scores=scores,
            poses=poses,
            metadata={
                "success": True,
                "error_message": None
            }
        )
        return result
    
    def create_task(self, ligand, receptor, configuration, metadata=None):
        """タスクを作成"""
        return DockingTask(
            id=str(uuid.uuid4()),
            ligand=ligand,
            receptor=receptor,
            configuration=configuration,
            metadata=metadata or {}
        )
    
    def validate_task(self, task):
        """タスクを検証"""
        return True
    
    def get_supported_parameters(self):
        """サポートされるパラメータ"""
        return ["exhaustiveness", "num_modes", "energy_range", "seed", "cpu"]
    
    def get_default_configuration(self):
        """デフォルト設定"""
        params = DockingParameters([
            DockingParameter("exhaustiveness", 8, ParameterType.INTEGER),
            DockingParameter("num_modes", 9, ParameterType.INTEGER),
            DockingParameter("energy_range", 3, ParameterType.FLOAT),
            DockingParameter("seed", 0, ParameterType.INTEGER),
            DockingParameter("cpu", 1, ParameterType.INTEGER)
        ])
        
        grid_box = GridBox(0.0, 0.0, 0.0, 20.0, 20.0, 20.0)
        
        return DockingConfiguration(
            grid_box=grid_box,
            parameters=params,
            name="Default Configuration",
            description="Default docking configuration"
        )
    
    def cancel_task(self, task):
        """タスクをキャンセル"""
        return True
    
    def get_result_for_task(self, task_id):
        """タスク結果を取得"""
        return None


def main():
    """メイン処理"""
    # 作業ディレクトリを設定
    work_dir = Path(".")
    input_dir = work_dir / "input" / "ALDR"
    output_dir = work_dir / "output" / "ALDR"
    os.makedirs(output_dir, exist_ok=True)
    
    # 入力ファイルのパス
    receptor_path = input_dir / "receptor.pdb"
    ligand_path = input_dir / "actives_final.sdf.gz"
    crystal_ligand_path = input_dir / "crystal_ligand.mol2"
    
    # ファイルの存在確認
    missing_files = []
    for path in [receptor_path, ligand_path, crystal_ligand_path]:
        if not path.exists():
            missing_files.append(path)
    
    if missing_files:
        print("以下のファイルが見つかりません：")
        for path in missing_files:
            print(f"  - {path}")
        print("必要なファイルを配置してから再度実行してください。")
        return
    
    print(f"タンパク質: {receptor_path}")
    print(f"リガンド: {ligand_path}")
    print(f"結晶リガンド: {crystal_ligand_path}")
    
    # タンパク質の準備
    protein = Protein(
        id="ALDR",
        path=str(receptor_path),
        structure=None,
        format=None
    )
    receptor = Receptor(protein=protein)
    print(f"タンパク質の準備が完了しました")
    
    # 1. 結晶リガンドファイル（MOL2）からグリッドボックスを計算
    print("\n方法1: 結晶リガンドからグリッドボックスを計算")
    try:
        # MOL2ファイルを解析して分子構造を作成
        crystal_structure = load_mol2_structure(str(crystal_ligand_path))
        print(f"結晶リガンドの原子数: {len(crystal_structure.atoms)}")
        
        # 計算された中心座標を表示
        center = crystal_structure.get_center_of_mass()
        print(f"計算されたドッキング中心: ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
        
        # 結晶リガンドからグリッドボックスを計算（パディング4Å）
        grid_box_ligand = GridBoxCalculator.calculate_from_structure(
            crystal_structure,
            padding=4.0,  # 周囲に4Åのパディングを追加
            default_box_size=None  # 実際のサイズから計算
        )
        
        # 計算されたボックスサイズを表示
        box_size = grid_box_ligand.get_size()
        print(f"計算されたボックスサイズ: {box_size[0]:.1f} x {box_size[1]:.1f} x {box_size[2]:.1f} Å")
        
        # 最小サイズの保証（10Å未満の場合は10Åに設定）
        min_size = 10.0
        if box_size[0] < min_size or box_size[1] < min_size or box_size[2] < min_size:
            print("ボックスサイズが小さすぎるため、最小サイズ（10Å x 10Å x 10Å）に調整します")
            size_x = max(box_size[0], min_size)
            size_y = max(box_size[1], min_size)
            size_z = max(box_size[2], min_size)
            grid_box_ligand = GridBox(
                center_x=grid_box_ligand.center_x,
                center_y=grid_box_ligand.center_y,
                center_z=grid_box_ligand.center_z,
                size_x=size_x,
                size_y=size_y,
                size_z=size_z
            )
    except Exception as e:
        print(f"結晶リガンドからのグリッドボックス計算中にエラーが発生しました: {e}")
        print("デフォルト値を使用します")
        grid_box_ligand = GridBox(0.0, 0.0, 0.0, 20.0, 20.0, 20.0)
    
    # 2. タンパク質構造からグリッドボックスを計算（代替方法）
    print("\n方法2: タンパク質全体からグリッドボックスを計算")
    grid_box_protein = None
    try:
        # 本来はここでタンパク質構造を解析しますが、今回は未実装のためスキップ
        # protein.calculate_grid_box() を使用する予定
        print("（この方法では構造解析未実装のため、デフォルト値を使用します）")
        grid_box_protein = GridBox(0.0, 0.0, 0.0, 20.0, 20.0, 20.0)
    except Exception as e:
        print(f"タンパク質からのグリッドボックス計算中にエラーが発生しました: {e}")
        grid_box_protein = GridBox(0.0, 0.0, 0.0, 20.0, 20.0, 20.0)
    
    # 結晶リガンドから計算したグリッドボックスを使用
    grid_box = grid_box_ligand
    
    # ドッキングパラメータの設定
    parameters = DockingParameters([
        DockingParameter("exhaustiveness", 8, ParameterType.INTEGER),
        DockingParameter("num_modes", 9, ParameterType.INTEGER),
        DockingParameter("energy_range", 3, ParameterType.FLOAT),
        DockingParameter("seed", 0, ParameterType.INTEGER),
        DockingParameter("cpu", 1, ParameterType.INTEGER)
    ])
    
    # ドッキング設定
    config = DockingConfiguration(
        grid_box=grid_box,
        parameters=parameters,
        name="ALDR Docking",
        description="ALDR タンパク質とリガンドのドッキング計算（自動計算グリッドボックス使用）"
    )
    
    # DockingServiceのインスタンス化
    docking_service = MockDockingService()
    
    # SDFファイルの読み込み
    try:
        suppl = Chem.SDMolSupplier(str(ligand_path))
        mols = [m for m in suppl if m is not None]
        print(f"\n読み込んだ化合物数: {len(mols)}")
        
        # 処理する化合物数の制限
        max_compounds = 5
        mols = mols[:min(max_compounds, len(mols))]
        
        # 各化合物に対してドッキングを実行
        for i, mol in enumerate(mols):
            mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Compound_{i+1}"
            print(f"\nリガンド {i+1}/{len(mols)}: {mol_name} を処理中...")
            
            # 化合物の準備
            compound = Compound(
                id=f"ligand_{i+1}",
                path=str(ligand_path),
                metadata={"name": mol_name}
            )
            # 名前を設定
            compound.name = mol_name
            
            # リガンドの準備
            ligand = Ligand(compound=compound, name=mol_name)
            ligand.set_prepared(True)
            
            # タスクの作成
            task = docking_service.create_task(ligand, receptor, config)
            
            # ドッキング実行
            result = docking_service.execute(task)
            
            # スコアの表示
            if result.scores:
                print(f"ドッキングスコア: {result.scores[0].value:.2f} kcal/mol")
            
            # 結果の保存（実際のPDBQT内容を書き込むにはポーズから文字列を生成する必要あり）
            result_path = output_dir / f"docking_result_{i+1}.pdbqt"
            try:
                with open(result_path, "w") as f:
                    # ポーズから内容を生成（今回はモックなので簡易的に）
                    pose_content = "REMARK Docking result\n"
                    pose_content += f"REMARK Score: {result.scores[0].value:.2f}\n"
                    # モックのPDBQT内容（実際は構造から正しく生成する必要あり）
                    for atom in result.poses[0].structure.atoms:
                        pose_content += f"ATOM  {'':>5} {'':2} {'':3} {'':>1} {'':>4}    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}\n"
                    f.write(pose_content)
                print(f"結果を {result_path} に保存しました")
            except Exception as e:
                print(f"結果ファイルの保存中にエラーが発生しました: {e}")
        
        print("\nすべてのドッキング計算が完了しました")
    
    except Exception as e:
        print(f"処理中にエラーが発生しました: {e}")


if __name__ == "__main__":
    main()