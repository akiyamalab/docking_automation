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
    
    chunk_size = 2  # チャンク分割の差異の、チャンクあたりの化合物数
    compound_set = CompoundSet(compound_path)
    compound_sets = compound_set.split_by_chunks(chunk_size)
    
    # ドッキング計算タスクの作成
    docking_tool = AutoDockVina()
    docking_tasks: list[Task] = []
    for i, split_compound_set in enumerate(compound_sets):
        # 各分割された化合物セットに対してタスクを作成
        task = Task.create(
            function=docking_tool.run_docking,
            args={
                "protein": Protein(protein_path),
                "compound_set": split_compound_set,
                "grid_box": GridBox(center=(15.0, 23.0, 36.0), size=(20.0, 20.0, 20.0)),
                "additional_params": AutoDockVinaParameters()
            },
            id=f"docking_task_{i}"
        )
        docking_tasks.append(task)

    # タスクの登録と並列計算の実行
    executor = DaskExecutor(scheduler_type="local")
    task_manager = TaskManager(executor=executor)
    for task in docking_tasks:
        task_manager.add_task(task)
    results = task_manager.execute_all()

    # 実行結果の統合
    all_docking_results = []
    for i, result in enumerate(results):
        # 結果を統合
        task_results = result.get_all()
        print(f"  - 結果数: {len(task_results)}")
        all_docking_results.extend(task_results)
    
    # 統合された結果からトップヒットを取得
    combined_results = DockingResultCollection()
    combined_results.extend(all_docking_results)
    
    # トップヒットを表示
    top_hits = combined_results.get_top(10)  # 上位10件を表示
    print(f"\n全タスクの結果から上位 {len(top_hits)} 件:")
    
    
    # 化合物セットのリストを作成
    compound_sets = [task.args["compound_set"] for task in docking_tasks]
    
    for j, hit in enumerate(top_hits):
        scores = hit.metadata["scores"]
        print(f"   全ポーズのスコア: {[float(s) for s in scores[0]]}")
    
    return results


if __name__ == "__main__":
    run_parallel_docking()