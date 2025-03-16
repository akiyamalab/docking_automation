#!/usr/bin/env python3
"""モックサービスを使用したバッチドッキングワークフローのテスト

このスクリプトは、モックサービスを使用してBatchDockingWorkflowの機能をテストします。
実際のドッキング計算は実行せず、シミュレーションのみを行います。
"""

import os
import sys
import uuid
import tempfile
import time
from pathlib import Path
from typing import List, Optional, Union, cast

# ライブラリのパスを追加
sys.path.insert(0, os.path.abspath('.'))

# モックシステムを使用（実際のドッキングツールは実行しない）
os.environ['MOCK_DOCKING'] = '1'

print("1. モジュールのインポート")
try:
    from docking_automation.application.workflows.batch_docking_workflow import BatchDockingWorkflow
    # dockingモジュールからのインポート（より簡潔に）
    from docking_automation.docking import (
        GridBox, DockingConfiguration, DockingParameters, DockingService,
        DockingTask, DockingResult, Ligand, Receptor, Score, ScoreType
    )
    # moleculeモジュールからのインポート（より簡潔に）
    from docking_automation.molecule import Compound, Protein, MoleculePreparationService
    print("モジュールのインポート成功")
except Exception as e:
    print(f"モジュールのインポートに失敗: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

# モック分子準備サービス
class MockMoleculePreparationService(MoleculePreparationService):
    """テスト用のモック分子準備サービス"""
    
    def prepare_ligand(self, compound, method=None):
        """リガンドを準備する（モック実装）"""
        print(f"リガンドを準備: {compound.path}")
        compound.add_metadata({"prepared": True})
        return compound
    
    def prepare_receptor(self, protein, method=None):
        """レセプターを準備する（モック実装）"""
        print(f"レセプターを準備: {protein.path}")
        protein.add_metadata({"prepared": True})
        return protein
    
    def calculate_properties(self, compound):
        """化合物の物理化学的特性を計算する（モック実装）"""
        return []
    
    def add_hydrogens(self, molecule, ph=7.4):
        """分子に水素原子を追加する（モック実装）"""
        return molecule
    
    def remove_water(self, protein):
        """タンパク質から水分子を除去する（モック実装）"""
        return protein
    
    def get_supported_methods(self):
        """サポートされている準備方法を返す（モック実装）"""
        return {"default": "モックメソッド"}

# モックドッキングサービス
class MockDockingService(DockingService):
    """テスト用のモックドッキングサービス"""
    
    def get_default_configuration(self):
        """デフォルトのドッキング設定を取得する（モック実装）"""
        return DockingConfiguration(
            grid_box=GridBox(
                center_x=0.0, center_y=0.0, center_z=0.0,
                size_x=20.0, size_y=20.0, size_z=20.0
            ),
            parameters=DockingParameters.from_dict({
                "exhaustiveness": 8,
                "num_modes": 9
            }),
            name="モックドッキング設定"
        )
    
    def validate_task(self, task):
        """タスクを検証する（モック実装）"""
        return True
    
    def get_supported_parameters(self):
        """サポートするパラメータを返す（モック実装）"""
        return {
            "exhaustiveness": "探索の徹底度（1-8）",
            "num_modes": "出力するポーズの数（1-10）",
            "energy_range": "エネルギー範囲（kcal/mol）",
            "seed": "乱数シード"
        }
    
    def get_tool_info(self):
        """ツール情報を取得する（モック実装）"""
        return {
            "name": "モックドッキングツール",
            "version": "1.0.0",
            "description": "テスト用のモックドッキングツール"
        }
    
    def cancel_task(self, task_id):
        """タスクをキャンセルする（モック実装）"""
        print(f"タスクをキャンセル: {task_id}")
        return True
    
    def get_result_for_task(self, task_id):
        """タスクの結果を取得する（モック実装）"""
        print(f"タスク結果を取得: {task_id}")
        return None  # 実際の結果は返さない（モック）
    
    def create_task(self, ligand, receptor, configuration, metadata=None):
        """ドッキングタスクを作成する（モック実装）"""
        return DockingTask(
            id=str(uuid.uuid4()),
            ligand=ligand,
            receptor=receptor,
            configuration=configuration,
            metadata=metadata or {}
        )
    
    def execute(self, task):
        """ドッキングタスクを実行する（モック実装）"""
        print(f"モックドッキングを実行: {task.id}")
        
        # モックスコア値を生成（リガンド名から決定的に値を生成）
        ligand_name = task.ligand.name
        # 文字列をハッシュ化して-5.0から-12.0の間のスコアを生成
        score_value = -5.0 - (hash(ligand_name) % 7000) / 1000.0
        
        # モック結果を生成
        return DockingResult(
            id=str(uuid.uuid4()),
            task=task,
            poses=[],  # 空のポーズリスト
            scores=[
                Score(
                    value=score_value,
                    score_type=ScoreType.BINDING_AFFINITY,
                    name="vina_score",
                    unit="kcal/mol"
                )
            ],
            created_at=time.time(),
            execution_time=0.5
        )

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
    
    # モックサービスの作成
    molecule_preparation_service = MockMoleculePreparationService()
    docking_service = MockDockingService()
    
    # ワークフローの作成（モックサービスを使用）
    workflow = BatchDockingWorkflow(
        molecule_preparation_service=molecule_preparation_service,
        docking_service=docking_service
    )
    
    # 入力パスをPath型に変換し、適切な型にキャスト
    receptor_path = Path(protein_file)
    ligand_paths = cast(List[Union[str, Path]], [Path(path) for path in ligand_files])
    
    # ワークフローを実行
    results = workflow.execute(  # execute()メソッドを使用
        receptor_path=receptor_path,
        ligand_paths=ligand_paths,  # すでにキャスト済み
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
        output_dir=Path(output_dir)
    )
    
    print(f"処理したリガンド数: {len(results)}")
    for i, result in enumerate(results[:3]):  # 最初の3つだけ表示
        score = result.get_best_score()
        score_value = score.value if score else "N/A"
        print(f"  結果 {i+1}: スコア={score_value}")
    
    # 出力ファイルの確認
    summary_file = os.path.join(output_dir, "batch_summary.txt")
    if os.path.exists(summary_file):
        print(f"\n出力サマリーファイル: {summary_file}")
        with open(summary_file, "r") as f:
            summary_content = f.read()
            print("サマリー内容（先頭3行）:")
            print("\n".join(summary_content.split("\n")[:3]) + "...")
    
except Exception as e:
    print(f"ワークフローの実行に失敗: {e}")
    import traceback
    traceback.print_exc()

print("\nテスト完了。")