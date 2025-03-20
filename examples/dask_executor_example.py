"""
DaskExecutorを使って複数のドッキング計算を並列実行するサンプル。

このサンプルでは、DaskExecutorを使って複数のDockingToolABC.run_docking()
メソッド呼び出しを並列に実行する方法を示します。
"""

from typing import List
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


def run_parallel_docking():
    """DaskExecutorを使って複数のドッキング計算を並列実行します。"""
    # 標準エラー出力をキャプチャするための設定
    # 元の標準エラー出力を保存
    original_stderr = sys.stderr
    # 標準エラー出力をStringIOにリダイレクト
    error_output = StringIO()
    sys.stderr = error_output
    
    try:
        # Initialize docking tool
        docking_tool = AutoDockVina()

        # Set up grid box and parameters
        grid_box = GridBox(center=(15.0, 23.0, 36.0), size=(20.0, 20.0, 20.0))
        additional_params = AutoDockVinaParameters(
            exhaustiveness=4,
            num_modes=3,
            seed=42
        )

        # 現在のスクリプトのディレクトリを取得
        script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
        
        # Load proteins and compound sets
        # タンパク質ファイルを使用する
        protein_path = script_dir / "input" / "ALDR" / "receptor.pdb"
        print(f"タンパク質ファイルの絶対パス: {protein_path.absolute()}")
        
        # 活性化合物を使用
        compound_path = script_dir / "input" / "ALDR" / "actives_final.sdf.gz"
        print(f"化合物ファイルの絶対パス: {compound_path.absolute()}")
        
        # Proteinオブジェクトの初期化
        try:
            print("タンパク質の前処理を開始します...")
            protein = Protein(protein_path)
            print("タンパク質の前処理が完了しました")
        except Exception as e:
            print(f"タンパク質の前処理中にエラーが発生しました: {e}")
            # エラー出力を即座に表示
            error_content = error_output.getvalue()
            if error_content:
                print("\n=== prepare_receptorのエラー出力 ===")
                print(error_content)
                print("=== エラー出力終了 ===\n")
            raise  # エラーを再度発生させて処理を中断
        
        # 化合物セットの初期化
        compound_set = CompoundSet(compound_path)
        
        # 化合物数を表示
        total_compounds = compound_set.get_compound_count()
        print(f"化合物数: {total_compounds}")
        
        # 処理する化合物の最大数を設定
        max_compounds = 3  # 3つの化合物のみを処理
        if max_compounds < total_compounds:
            print(f"最初の{max_compounds}件の化合物のみを処理します")
        
        # Create tasks for each compound
        docking_tasks = []
        for i in range(max_compounds):
            # 各化合物に対してタスクを作成
            task = Task.create(
                function=docking_tool.run_docking,
                args={
                    "protein": protein,
                    "compound_set": compound_set,
                    "grid_box": grid_box,
                    "additional_params": additional_params
                },
                id=f"docking_task_{i}"
            )
            docking_tasks.append(task)
    except Exception as e:
        print(f"エラーが発生しました: {e}")
        # エラー出力を表示
        error_content = error_output.getvalue()
        if error_content:
            print("\n=== 標準エラー出力 ===")
            print(error_content)
            print("=== エラー出力終了 ===\n")
        raise
    finally:
        # 標準エラー出力を元に戻す
        sys.stderr = original_stderr

    # Choose executor based on environment
    if os.environ.get("SLURM_JOB_ID"):
        # Slurm environment
        executor = DaskExecutor(
            scheduler_type="slurm",
            jobs=10,
            cores=4,
            memory="16GB",
            walltime="02:00:00",
            queue="normal"
        )
        print("Using Slurm environment")
    else:
        # Local environment
        executor = DaskExecutor(scheduler_type="local")
        print("Using local environment")

    # Set up task manager and add tasks
    task_manager = TaskManager(executor=executor)
    for task in docking_tasks:
        task_manager.add_task(task)

    # Execute all tasks in parallel
    print(f"Executing {len(docking_tasks)} docking tasks in parallel...")
    results = task_manager.execute_all()

    # Process results
    for i, result in enumerate(results):
        if result is not None:
            print(f"Task {i} result: {result}")
            # 結果の詳細を表示
            if hasattr(result, 'get_top'):
                top_hits = result.get_top(3)  # 上位3件を表示
                print(f"\nTop {len(top_hits)} hits:")
                for j, hit in enumerate(top_hits):
                    print(f"{j+1}. Score: {hit.docking_score}, Compound: {hit.compound_set_id}_{hit.compound_index}")
                    print(f"   Pose file: {hit.result_path}")
        else:
            print(f"Task {i} result: None")

    return results


if __name__ == "__main__":
    run_parallel_docking()