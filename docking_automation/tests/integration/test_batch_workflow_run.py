#!/usr/bin/env python3
"""バッチドッキングワークフローの実行テスト

このスクリプトは、BatchDockingWorkflowが複数のリガンドファイルを
処理できることを確認するためのテストを実行します。
"""

import os
import sys
import tempfile
import uuid
import time
from pathlib import Path
from typing import List, Union, Optional, cast

# ライブラリのパスを追加
sys.path.insert(0, os.path.abspath('.'))

# モックシステムを使用（実際のドッキングツールは実行しない）
os.environ['MOCK_DOCKING'] = '1'

print("1. モジュールのインポート")
try:
    from docking_automation.application.workflows.batch_docking_workflow import BatchDockingWorkflow
    # dockingモジュールからのインポート（より簡潔に）
    from docking_automation.docking import (
        GridBox, DockingService, Ligand, Receptor, DockingTask, DockingResult,
        Score, ScoreType, DockingConfiguration, DockingParameters
    )
    # moleculeモジュールからのインポート（より簡潔に）
    from docking_automation.molecule import Compound, Protein, MoleculePreparationService
    print("モジュールのインポート成功")
except Exception as e:
    print(f"モジュールのインポートに失敗: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# モックサービスの定義
class MockMoleculePreparationService(MoleculePreparationService):
    """モック分子準備サービス（テスト用）"""
    def prepare_ligand(self, compound, method=None):
        print(f"モック：リガンドを準備: {compound.path}")
        return compound
        
    def prepare_receptor(self, protein, method=None):
        print(f"モック：レセプターを準備: {protein.path}")
        return protein
        
    def calculate_properties(self, compound):
        return []
        
    def add_hydrogens(self, molecule, ph=7.4):
        return molecule
        
    def remove_water(self, protein):
        return protein
        
    def get_supported_methods(self):
        return {"default": "モックメソッド"}

class MockDockingService(DockingService):
    """モックドッキングサービス（テスト用）"""
    def get_default_configuration(self):
        return DockingConfiguration(
            grid_box=GridBox(center_x=0.0, center_y=0.0, center_z=0.0, size_x=20.0, size_y=20.0, size_z=20.0),
            parameters=DockingParameters.from_dict({"exhaustiveness": 8}),
            name="モック設定"
        )
        
    def create_task(self, ligand, receptor, configuration, metadata=None):
        return DockingTask(
            id=str(uuid.uuid4()),
            ligand=ligand,
            receptor=receptor,
            configuration=configuration,
            metadata=metadata or {}
        )
        
    def execute(self, task):
        print(f"モック：ドッキングを実行: {task.id}")
        score_value = -8.5  # モックスコア値
        return DockingResult(
            id=str(uuid.uuid4()),
            task=task,
            poses=[],
            scores=[Score(value=score_value, score_type=ScoreType.BINDING_AFFINITY, name="mock_score", unit="kcal/mol")],
            created_at=time.time(),
            execution_time=0.1
        )
        
    def validate_task(self, task):
        return True
        
    def get_supported_parameters(self):
        return {}
        
    def get_tool_info(self):
        return {"name": "MockDocking"}
        
    def cancel_task(self, task_id):
        return True
        
    def get_result_for_task(self, task_id):
        return None

print("\n2. テスト準備")
# テストに使用するファイルパスを設定
input_dir = "input_ligands"
if os.path.exists(input_dir):
    ligand_files = [
        os.path.join(input_dir, file)
        for file in os.listdir(input_dir)
        if file.endswith((".sdf", ".pdb"))
    ][:3]  # 最初の3つだけ使用
    
    print(f"テスト用リガンドファイル: {ligand_files}")
else:
    # 入力ディレクトリがなければテスト用ファイルを作成
    with tempfile.TemporaryDirectory() as temp_dir:
        ligand_files = [os.path.join(temp_dir, f"test_ligand_{i}.sdf") for i in range(3)]
        for file in ligand_files:
            with open(file, 'w') as f:
                f.write("MOCK SDF CONTENT\n")
        
        print(f"テスト用一時ファイル作成: {ligand_files}")

# タンパク質ファイルの設定
protein_file = "protein.pdbqt"
if not os.path.exists(protein_file):
    # テスト用のモックファイルを作成
    with open(protein_file, 'w') as f:
        f.write("MOCK PDBQT CONTENT\n")
    print(f"テスト用モックタンパク質ファイル作成: {protein_file}")
    
print("\n3. ワークフローの実行")
try:
    # 出力ディレクトリの設定
    output_dir = "batch_test_results"
    os.makedirs(output_dir, exist_ok=True)
    
    # サービスを作成
    molecule_preparation_service = MockMoleculePreparationService()
    docking_service = MockDockingService()
    
    # ワークフローの作成（必要なサービスを渡す）
    workflow = BatchDockingWorkflow(
        molecule_preparation_service=molecule_preparation_service,
        docking_service=docking_service
    )
    
    # Path型に変換
    receptor_path = Path(protein_file)
    # 正確な型にキャスト（List[Union[str, Path]]）
    ligand_paths = cast(List[Union[str, Path]], [Path(p) for p in ligand_files])
    output_path = Path(output_dir)
    
    # ワークフローを実行
    results = workflow.execute(
        receptor_path=receptor_path,
        ligand_paths=ligand_paths,
        config={
            "grid_box": {
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 20.0,
                "size_y": 20.0,
                "size_z": 20.0
            }
        },
        output_dir=output_path
    )
    
    print(f"処理したリガンド数: {len(results)}")
    for i, result in enumerate(results[:3]):  # 最初の3つだけ表示
        score = result.get_best_score()
        score_value = score.value if score else "N/A"
        print(f"  結果 {i+1}: スコア={score_value}")
    
except Exception as e:
    print(f"ワークフローの実行に失敗: {e}")
    import traceback
    traceback.print_exc()

print("\nテスト完了。")