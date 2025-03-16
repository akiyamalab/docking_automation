#!/usr/bin/env python3
"""バッチドッキングワークフローと結果解析のテスト

このスクリプトは、修正したBatchDockingWorkflowの動作と
ResultAnalysisServiceのヒストグラム生成機能をテストします。
"""

import os
import sys
import uuid
import tempfile
from pathlib import Path
import matplotlib
# GUIを使用しない（サーバー環境用）
matplotlib.use('Agg')

# ライブラリのパスを追加
sys.path.insert(0, os.path.abspath('.'))

# 必要なクラスをインポート
print("1. モジュールのインポート")
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.entity.docking_result import DockingResult
from docking_automation.docking.value_object.score import Score, ScoreType
from docking_automation.docking.value_object.docking_parameter import DockingParameters
from docking_automation.application.services.result_analysis_service import ResultAnalysisService

print("モジュールのインポート成功")

# テスト用のモックデータを作成
def create_mock_result(compound_id, score_value):
    """モックのドッキング結果を作成"""
    from docking_automation.docking.entity.docking_task import DockingTask
    
    # 化合物を作成
    compound = Compound(id=compound_id, path=f"mock/compound_{compound_id}.sdf")
    
    # タンパク質を作成
    protein = Protein(id="P001", path="mock/protein.pdb")
    
    # リガンドとレセプターを作成
    ligand = Ligand(compound=compound)
    receptor = Receptor(protein=protein)
    
    # パラメータとタスクを作成
    parameters = DockingParameters.from_dict({
        'exhaustiveness': 8,
        'num_modes': 9,
    })
    
    task = DockingTask(
        id=f"T{compound_id}",
        ligand=ligand,
        receptor=receptor,
        configuration=DockingConfiguration(
            grid_box=GridBox(center_x=0, center_y=0, center_z=0, size_x=10, size_y=10, size_z=10),
            parameters=parameters,
            name="mock_config"
        )
    )
    
    # スコアを作成
    score = Score(
        value=score_value,
        name="vina_score",
        score_type=ScoreType.BINDING_AFFINITY,
        unit="kcal/mol"
    )
    
    # 結果を作成（必須パラメータを追加）
    import time
    result = DockingResult(
        id=f"R{compound_id}",
        task=task,
        poses=[],  # 空のポーズリスト
        scores=[score],
        created_at=time.time(),
        execution_time=1.0
    )
    
    return result

# 複数の結果を作成してテスト
print("\n2. 結果解析サービスのテスト")
try:
    # 一時ディレクトリを作成
    with tempfile.TemporaryDirectory() as temp_dir:
        # モックの結果を10個作成
        mock_results = [
            create_mock_result(f"C{i:03d}", -8.0 - i * 0.5) for i in range(10)
        ]
        
        # 結果解析サービスを作成
        analysis_service = ResultAnalysisService()
        
        # スコア統計を計算
        scores = [r.scores[0].value for r in mock_results]
        stats = analysis_service.calculate_score_statistics(scores)
        print(f"スコア統計: {stats}")
        
        # ヒストグラムを生成
        histogram_path = os.path.join(temp_dir, "score_histogram.png")
        analysis_service.generate_score_histogram(
            results=mock_results,
            output_path=histogram_path
        )
        
        print(f"ヒストグラム生成: {os.path.exists(histogram_path)}")
        
        # ランキングのテスト
        ranked_results = analysis_service.rank_results_by_score(mock_results)
        print(f"ランキング結果（上位3件）:")
        for i, result in enumerate(ranked_results[:3]):
            print(f"  {i+1}. Score: {result.scores[0].value}")
        
        # 比較のテスト
        mock_results2 = [
            create_mock_result(f"D{i:03d}", -7.0 - i * 0.3) for i in range(10)
        ]
        
        comparison_path = os.path.join(temp_dir, "comparison.png")
        comparison = analysis_service.compare_results(
            mock_results, mock_results2,
            output_path=comparison_path
        )
        
        print(f"比較結果: {comparison['mean_diff']:.2f}")
        print(f"比較グラフ生成: {os.path.exists(comparison_path)}")
        
except Exception as e:
    print(f"結果解析サービスのテストに失敗: {e}")
    import traceback
    traceback.print_exc()

print("\nテスト完了。")