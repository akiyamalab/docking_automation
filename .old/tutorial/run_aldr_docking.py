#!/usr/bin/env python
# coding: utf-8

"""
ALDRタンパク質のドッキング計算を実行するスクリプト
crystal_ligand.mol2を参照してドッキング中心とボックスサイズを設定
"""

import os
import sys
from pathlib import Path
import random
import uuid
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
from docking_automation.docking.value_object.score import ScoreType
from docking_automation.docking.value_object.pose import Pose
from docking_automation.docking.value_object.score import Score
from docking_automation.docking.service.docking_service import DockingService

# モック実装を作成（本来はVinaのサービスを使用するが、ここではモックで代用）
class MockDockingService(DockingService):
    """デモ用のモックドッキングサービス"""
    def execute(self, task: DockingTask) -> DockingResult:
        """ドッキング計算を実行（モック）"""
        print(f"[モック] {task.ligand.name}のドッキング計算を実行中...")
        
        # ランダムなスコアを生成
        from docking_automation.docking.value_object.score import ScoreType
        scores = [Score(
            value=random.uniform(-12.0, -5.0),
            score_type=ScoreType.BINDING_AFFINITY,
            name=f"score_{i+1}"
        ) for i in range(3)]
        
        # モックのポーズを生成
        poses = []
        for i in range(3):
            # MoleculeStructureとAtomを作成してポーズにセット
            from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom
            
            # モック用の原子を作成
            atoms = [
                Atom(
                    atom_id=j+1,
                    element="C",
                    x=task.configuration.grid_box.center_x + (j * 0.5),
                    y=task.configuration.grid_box.center_y + (j * 0.2),
                    z=task.configuration.grid_box.center_z + (j * 0.3),
                )
                for j in range(10)
            ]
            
            structure = MoleculeStructure(atoms=atoms, bonds=[])
            poses.append(Pose(structure=structure, rank=i+1))
        
        # 結果を返す（created_atは浮動小数点のタイムスタンプが必要）
        import time
        return DockingResult(
            id=f"result_{task.id}",
            task=task,
            scores=scores,
            poses=poses,
            created_at=time.time()
        )
    
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
        
        return DockingConfiguration(grid_box=grid_box, parameters=params)
    
    def cancel_task(self, task):
        """タスクをキャンセル"""
        return True
    
    def get_result_for_task(self, task_id):
        """タスク結果を取得"""
        return None

# 作業ディレクトリを設定
work_dir = Path(".")
input_dir = work_dir / "input" / "ALDR"
output_dir = work_dir / "output" / "ALDR"
os.makedirs(output_dir, exist_ok=True)

# 入力ファイルのパス
receptor_path = input_dir / "receptor.pdb"
ligand_path = input_dir / "actives_final.sdf.gz"

# ドッキングの設定
docking_center = (16.645, -6.714, 14.124)  # crystal_ligand.mol2から計算
box_size = (20.0, 20.0, 20.0)  # 各方向に20Åのボックスサイズ

print(f"タンパク質: {receptor_path}")
print(f"リガンド: {ligand_path}")
print(f"ドッキング中心: {docking_center}")
print(f"ボックスサイズ: {box_size}")

# タンパク質の準備
protein = Protein(
    id="ALDR",
    path=str(receptor_path),
    structure=None,
    format=None
)
receptor = Receptor(protein=protein)
print(f"タンパク質の準備が完了しました")

# グリッドボックスの設定
grid_box = GridBox(
    center_x=docking_center[0],
    center_y=docking_center[1],
    center_z=docking_center[2],
    size_x=box_size[0],
    size_y=box_size[1],
    size_z=box_size[2]
)

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
    description="ALDR タンパク質とリガンドのドッキング計算"
)

# DockingServiceのインスタンス化
docking_service = MockDockingService()

# 圧縮されたSDFファイルを展開して読み込む
import gzip
import tempfile

# 一時ファイルに展開
temp_sdf = os.path.join(output_dir, "actives_final.sdf")
with gzip.open(ligand_path, 'rb') as f_in:
    with open(temp_sdf, 'wb') as f_out:
        f_out.write(f_in.read())

# SDFファイルの読み込み
suppl = Chem.SDMolSupplier(temp_sdf)
mols = [m for m in suppl if m is not None]
print(f"読み込んだ化合物数: {len(mols)}")

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
        structure=mol,
        metadata={"name": mol_name}
    )
    
    # リガンドの準備
    ligand = Ligand(compound=compound, name=mol_name)
    ligand.set_prepared(True)
    
    # タスクの作成
    task = docking_service.create_task(ligand, receptor, config)
    
    # ドッキング実行
    result = docking_service.execute(task)
    
    print(f"ドッキングスコア: {result.scores[0].value:.2f} kcal/mol")
    
    # 結果の保存
    result_path = output_dir / f"docking_result_{i+1}.pdbqt"
    with open(result_path, "w") as f:
        # MoleculeStructureからPDBQT形式のテキストを生成
        structure = result.poses[0].structure
        atoms = structure.atoms
        f.write("MODEL 1\n")
        for atom in atoms:
            pdbqt_line = f"ATOM  {atom.atom_id:5d}  {atom.element:<3s} LIG A   1    {atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}  1.00  0.00    0.000 {atom.element}\n"
            f.write(pdbqt_line)
        f.write("ENDMDL\n")
    
    print(f"結果を {result_path} に保存しました")

print("\nすべてのドッキング計算が完了しました")