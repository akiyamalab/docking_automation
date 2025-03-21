"""
DaskExecutorを使って複数のドッキング計算を並列実行するサンプル。

このサンプルでは、DaskExecutorを使って複数のDockingToolABC.run_docking()
メソッド呼び出しを並列に実行する方法を示します。

化合物セットを複数のグループに分割し、それぞれを1つのタスクとして並列処理します。
"""

from typing import List, Tuple
import os
import sys
from io import StringIO
from pathlib import Path

# OpenBabelのwarningを完全に抑制するための設定
os.environ["BABEL_QUIET"] = "1"

from docking_automation.infrastructure.executor.task import Task
from docking_automation.infrastructure.executor.task_manager import TaskManager
from docking_automation.infrastructure.executor.dask_executor import DaskExecutor
from docking_automation.docking.autodockvina_docking import AutoDockVina, AutoDockVinaParameters
from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.protein import Protein
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.docking.docking_result_collection import DockingResultCollection

# 各種ファイルのパスをハードコーディング
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
protein_path = script_dir / "input" / "ALDR" / "receptor.pdb"
compound_path = script_dir / "input" / "ALDR" / "actives_subset.sdf"

def run_parallel_docking():
    """
    DaskExecutorとTaskManagerを使って複数のドッキング計算を並列実行します。
    
    このメソッドでは、以下の手順で並列処理を行います：
    1. 化合物セットを複数のチャンクに分割
    2. 各チャンクに対してドッキングタスクを作成
    3. DaskExecutorとTaskManagerを使用して並列実行
    4. 結果を統合して上位ヒットを表示
    
    並列実行により、複数のCPUコアを活用して処理時間を短縮します。
    """
    
    # Set up grid box and parameters
    grid_box = GridBox(center=(15.0, 23.0, 36.0), size=(20.0, 20.0, 20.0))
    additional_params = AutoDockVinaParameters()

    # 現在のスクリプトのディレクトリを取得
    print(f"タンパク質ファイルの絶対パス: {protein_path.absolute()}")
    print(f"化合物ファイルの絶対パス: {compound_path.absolute()}")
    
    
    protein = Protein(protein_path)
    compound_set = CompoundSet(compound_path)
    
    chunk_size = 2  # チャンク分割の差異の、チャンクあたりの化合物数
    compound_sets = compound_set.split_by_chunks(chunk_size)
    
    docking_tool = AutoDockVina()
    docking_tasks: list[Task] = []
    for i, split_compound_set in enumerate(compound_sets):
        # 各分割された化合物セットに対してタスクを作成
        task = Task.create(
            function=docking_tool.run_docking,
            args={
                "protein": protein,
                "compound_set": split_compound_set,
                "grid_box": grid_box,
                "additional_params": additional_params
            },
            id=f"docking_task_{i}"
        )
        docking_tasks.append(task)

    # DaskExecutorの初期化
    executor = DaskExecutor(scheduler_type="local")
    
    # TaskManagerの初期化
    task_manager = TaskManager(executor=executor)
    
    # タスクをTaskManagerに追加
    for task in docking_tasks:
        task_manager.add_task(task)
    
    # タスクの並列実行
    print(f"並列計算を開始します（タスク数: {len(docking_tasks)}）...")
    results = task_manager.execute_all()

    # Process results
    all_docking_results = []
    for i, result in enumerate(results):
        # DockingResultCollectionの場合、結果を統合
        # TODO: get_all と get_all_results が両方あるのは不適切なので、統一する
        if hasattr(result, 'get_all'):
            task_results = result.get_all()
            print(f"  - 結果数: {len(task_results)}")
            all_docking_results.extend(task_results)
        elif hasattr(result, 'get_all_results'):
            task_results = result.get_all_results()
            print(f"  - 結果数: {len(task_results)}")
            all_docking_results.extend(task_results)
        else:
            print(f"  - 予期しない結果タイプ: {type(result)}")
    
    # 統合された結果からトップヒットを取得
    combined_results = DockingResultCollection()
    combined_results.extend(all_docking_results)
    
    # トップヒットを表示
    top_hits = combined_results.get_top(10)  # 上位10件を表示
    print(f"\n全タスクの結果から上位 {len(top_hits)} 件:")
    for j, hit in enumerate(top_hits):
        # TODO: このあたりの処理はある程度の部分について DockingResult に適切に実装されているべき
        
        # 化合物IDを抽出（actives_subset_range_X_Yから元のインデックスを取得）
        compound_id = hit.compound_set_id
        original_index = hit.compound_index
        
        # 元のインデックスを取得するための処理
        # 各タスクの実行結果を確認
        for i, task_result in enumerate(results):
            if task_result is not None:
                
                # タスクの化合物セットのIDを取得
                # 対応するタスクを探す
                for task in docking_tasks:
                    task_compound_set = task.args["compound_set"]
                    task_compound_set_id = task_compound_set.id
                    # 注：この部分は上のループ内に移動しました
                
                # タスクの化合物セットのIDと一致する場合
                if task_compound_set_id == compound_id:
                    # インデックス範囲を取得
                    properties = task_compound_set.get_properties()
                    if "index_range" in properties:
                        index_range = properties["index_range"]
                        # 元のインデックスを計算
                        original_index = index_range["start"] + hit.compound_index
                        break
        
        # 結果を表示
        print(f"{j+1}. Score: {hit.docking_score}, Compound: actives_subset_{original_index}")
        
        # 化合物の詳細情報を表示
        # これも DockingResult の責務にすべき
        if hit.metadata and "compound_name" in hit.metadata:
            print(f"   化合物名: {hit.metadata['compound_name']}")
        else:
            # メタデータがない場合は化合物IDを表示
            print(f"   化合物ID: actives_subset_{original_index}")
            
        # 必要に応じてメタデータの他の情報も表示
        if hit.metadata:
            if "scores" in hit.metadata and len(hit.metadata["scores"]) > 0:
                # 全てのポーズのスコアを表示
                scores = hit.metadata["scores"]
                print(f"   全ポーズのスコア: {[float(s) for s in scores[0]]}")
    
    return results


if __name__ == "__main__":
    run_parallel_docking()