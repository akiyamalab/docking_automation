"""
DaskExecutorを使って複数のドッキング計算を並列実行するサンプル。

このサンプルでは、DaskExecutorを使って複数のDockingToolABC.run_docking()
メソッド呼び出しを並列に実行する方法を示します。

Alphafoldで作成されたタンパク質構造をセグメンテーションし、
それぞれの構造に対してドッキング計算を並列実行します。
"""

import os
import shutil
from pathlib import Path
from typing import Any, Union

import numpy as np

# OpenBabelのwarningを完全に抑制するための設定
os.environ["BABEL_QUIET"] = "1"

from docking_automation.docking import (
    AutoDockVina,
    AutoDockVinaParameters,
    DockingResultCollection,
    GridBox,
)
from docking_automation.domain.services.protein_segmentation_service import (
    ProteinSegmentationService,
)
from docking_automation.infrastructure.executor import DaskExecutor, Task, TaskManager
from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
    DockingResultRepositoryFactory,
    RepositoryType,
)
from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
    HDF5DockingResultRepository,
)

# 時間計測用デコレータをインポート
from docking_automation.infrastructure.utilities.time_utils import (
    measure_execution_time,
)
from docking_automation.molecule import CompoundSet, Protein

# 各種ファイルのパスをハードコーディング
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
alphafold_protein_paths = [
    script_dir / "input" / "alphafold" / "AF-A0A0A2JW93-F1-model_v4.pdb",
    script_dir / "input" / "alphafold" / "AF-A0A0H2UNG0-F1-model_v4.pdb",
]
compound_path = script_dir / "input" / "ALDR" / "actives_subset.sdf"
output_dir = script_dir / "output" / "alphafold_segmentation"
repository_dir = output_dir / "repository"  # リポジトリディレクトリを追加
sdf_output_dir = output_dir / "sdf_exports"  # SDFエクスポート用ディレクトリを追加
protein_output_dir = sdf_output_dir / "proteins"  # タンパク質構造ファイル用ディレクトリを追加

# 出力ディレクトリとリポジトリディレクトリの作成
os.makedirs(output_dir, exist_ok=True)
os.makedirs(repository_dir, exist_ok=True)
os.makedirs(sdf_output_dir, exist_ok=True)
os.makedirs(protein_output_dir, exist_ok=True)


def _segment_protein(protein_path: Path, output_dir: Path) -> list[Protein]:
    """
    タンパク質構造をセグメンテーションする

    Args:
        protein_path: タンパク質構造ファイルのパス
        output_dir: 出力ディレクトリのパス

    Returns:
        セグメンテーションされたタンパク質構造のリスト
    """
    # 出力ディレクトリの作成
    os.makedirs(output_dir, exist_ok=True)

    # Proteinオブジェクトの作成
    protein = Protein.create(path=protein_path)
    print(str(protein))

    # ProteinSegmentationServiceの作成と実行
    # 進捗状況を標準出力に表示するコールバック関数
    def progress_callback(message: str):
        print(message)

    service = ProteinSegmentationService(progress_callback=progress_callback)
    segmented_proteins = service.segment(protein, None, output_dir)

    # 結果の表示
    print(service.get_segmentation_summary(segmented_proteins, protein))

    return segmented_proteins


def copy_protein_structure(protein: "Protein", output_dir: Path, name_prefix: str = "") -> Path:
    """
    タンパク質構造ファイルを指定されたディレクトリにコピーする。

    Args:
        protein: コピーするタンパク質
        output_dir: 出力先ディレクトリ
        name_prefix: ファイル名の接頭辞（オプション）

    Returns:
        コピー先のファイルパス
    """
    # 出力ディレクトリが存在しない場合は作成
    os.makedirs(output_dir, exist_ok=True)

    # 出力ファイル名を生成
    if name_prefix:
        output_filename = f"{name_prefix}_{protein.path.name}"
    else:
        output_filename = protein.path.name

    # 出力先のパスを生成
    output_path = output_dir / output_filename

    # ファイルをコピー
    shutil.copy2(protein.path, output_path)
    print(f"[Protein:{protein.id}] 構造ファイルをコピーしました: {output_path}")

    return output_path


@measure_execution_time
def run_docking(
    protein: "Protein",
    compound_set: Union["CompoundSet"],
    grid_box: "GridBox",
    additional_params: Any,
    repository_dir: Path = None,
) -> "DockingResultCollection":
    """
    ドッキング計算を実行する。計算結果の再利用機能を使用する。

    Args:
        protein: タンパク質
        compound_set: 化合物セット
        grid_box: グリッドボックス
        additional_params: 追加パラメータ
        repository_dir: リポジトリディレクトリのパス（指定しない場合は再利用しない）

    Returns:
        ドッキング結果のコレクション
    """
    # CompoundSetのインデックス範囲を取得
    start_index = 0
    try:
        properties = compound_set.get_properties()
        index_range = properties.get("index_range")
        if index_range is not None:
            start_index = index_range["start"]
    except Exception as e:
        print(f"[Protein:{protein.id}] インデックス範囲の取得中にエラーが発生しました: {e}")

    # リポジトリの初期化（指定されている場合）
    repository = None
    if repository_dir is not None:
        try:
            # HDF5リポジトリを初期化
            repository_factory = DockingResultRepositoryFactory()
            repository = repository_factory.create(
                repository_type=RepositoryType.HDF5,
                base_directory=repository_dir,
                config={"mode": "append"}
            )
            print(f"[Protein:{protein.id}] リポジトリを初期化しました: {repository_dir}")
        except Exception as e:
            print(f"[Protein:{protein.id}] リポジトリの初期化中にエラーが発生しました: {e}")
            # リポジトリの初期化に失敗した場合は、再利用せずに続行
            repository = None

    # ドッキング計算を実行
    docking_tool = AutoDockVina()
    
    # リポジトリが指定されている場合は、再利用機能を使用
    if repository is not None and isinstance(repository, HDF5DockingResultRepository):
        print(f"[Protein:{protein.id}] 計算結果の再利用機能を使用します")
        results = docking_tool.run_docking_with_reuse(
            protein=protein,
            compound_set=compound_set,
            grid_box=grid_box,
            additional_params=additional_params,
            repository=repository,
        )
    else:
        # 通常のドッキング計算を実行
        results = docking_tool.run_docking(
            protein=protein,
            compound_set=compound_set,
            grid_box=grid_box,
            additional_params=additional_params,
        )

    # compound_indexを修正
    for i, result in enumerate(results):
        # DockingResultクラスのcompound_indexフィールドを直接修正
        result.compound_index = start_index + i
        
        # 再利用されたかどうかを示すメタデータを追加
        if "reused" in result.metadata and result.metadata["reused"]:
            print(f"[Protein:{protein.id}] 化合物 {result.compound_index} の結果を再利用しました（スコア: {result.docking_score}）")

    # NumPy配列をリストに変換
    for result in results:
        for key, value in result.metadata.items():
            if isinstance(value, np.ndarray):
                result.metadata[key] = value.tolist()
            elif isinstance(value, list) and len(value) > 0 and isinstance(value[0], np.ndarray):
                result.metadata[key] = [arr.tolist() for arr in value]

    return results


def run_parallel_docking():
    """
    DaskExecutorとTaskManagerを使って複数のドッキング計算を並列実行します。

    このメソッドでは、以下の手順で並列処理を行います：
    1. Alphafoldタンパク質構造をセグメンテーション
    2. 化合物セットを複数のチャンクに分割
    3. 各セグメントと各チャンクの組み合わせに対してドッキングタスクを作成
    4. DaskExecutorとTaskManagerを使用して並列実行
    5. 結果を統合して上位ヒットを表示

    並列実行により、複数のCPUコアを活用して処理時間を短縮します。

    Returns:
        list: 各タスクの実行結果
    """
    # タンパク質構造のセグメンテーション
    all_segmented_proteins = []
    protein_names = []
    for protein_path in alphafold_protein_paths:
        protein_name = protein_path.stem
        protein_output_dir = output_dir / protein_name
        print(f"\n=== タンパク質 {protein_name} の処理を開始 ===")

        segmented_proteins = _segment_protein(protein_path, protein_output_dir)
        all_segmented_proteins.extend(segmented_proteins)
        protein_names.extend([protein_name] * len(segmented_proteins))

    # 化合物セットの読み込みと分割
    chunk_size = 1  # チャンク分割の、チャンクあたりの化合物数
    compound_set = CompoundSet(compound_path)
    compound_sets: list[CompoundSet] = compound_set.split_by_chunks(chunk_size)

    # ドッキング計算タスクの作成
    docking_tasks: list[Task] = []

    # 各セグメントに対してドッキング範囲の推定とドッキングタスクの作成
    grid_boxes = {}
    print("\n=== タンパク質・セグメントごとのドッキング範囲推定 ===")

    for i, (protein, protein_name) in enumerate(zip(all_segmented_proteins, protein_names)):
        # GridBoxの取得
        print(f"\nタンパク質: {protein_name}, セグメント {i + 1} ({protein.id}) のドッキング範囲を推定中...")
        grid_box = GridBox.from_fpocket(protein)

        # グリッドボックス情報の表示
        print(str(grid_box))

        # グリッドボックスを保存
        grid_boxes[(protein_name, i)] = grid_box

        for j, split_compound_set in enumerate(compound_sets):
            # 各分割された化合物セットに対してタスクを作成
            task = Task.create(
                function=run_docking,  # 関数名を変更
                args={
                    "protein": protein,
                    "compound_set": split_compound_set,
                    "grid_box": grid_box,
                    "additional_params": AutoDockVinaParameters(),
                    "repository_dir": repository_dir,  # リポジトリディレクトリを追加
                },
                id=f"docking_task_protein_{protein_name}_segment_{i}_chunk_{j}",
            )
            docking_tasks.append(task)

    # タスクの登録と並列計算の実行
    executor = DaskExecutor(scheduler_type="local", n_workers=None)
    task_manager = TaskManager(executor=executor)
    for task in docking_tasks:
        task_manager.add_task(task)

    # リポジトリ情報をタスクマネージャに渡す（スケジューラ側での一元的な保存のため）
    results = task_manager.execute_all(
        {
            "repository_type": RepositoryType.HDF5,
            "base_directory": repository_dir,
            "config": {"mode": "append"},
        }
    )

    # 各タンパク質・セグメントごとの結果を保持
    protein_segment_results = {}
    for i, task in enumerate(docking_tasks):
        task_id = task.id
        protein = task.args["protein"]

        # タスクIDからタンパク質名とセグメント情報を抽出
        protein_name = task_id.split("protein_")[1].split("_segment_")[0]
        segment_index = int(task_id.split("_segment_")[1].split("_chunk_")[0])

        # 結果を格納するキー
        result_key = (protein_name, segment_index)

        if result_key not in protein_segment_results:
            protein_segment_results[result_key] = {
                "protein": protein,
                "protein_name": protein_name,
                "segment_index": segment_index,
                "results": DockingResultCollection(),
            }

        # 結果を対応するタンパク質・セグメントに追加
        task_results = results[i].get_all()

        # 各ドッキング結果にタスクIDを追加
        for result in task_results:
            result.metadata["task_id"] = task_id

        protein_segment_results[result_key]["results"].extend(task_results)

    # 各タンパク質・セグメントごとのトップヒットを表示
    for (protein_name, segment_index), data in protein_segment_results.items():
        protein = data["protein"]
        grid_box = grid_boxes.get((protein_name, segment_index))

        print(f"\nタンパク質: {protein_name}, セグメント {segment_index + 1} ({protein.id}) の結果:")
        if grid_box:
            print(f"  ドッキング範囲: {str(grid_box)}")

        # 上位ヒットの表示
        print(data["results"].format_top_hits(5))

        # タンパク質・セグメントごとの結果をSDFファイルにエクスポート
        try:

            # SDFファイルにエクスポート
            protein_segment_id = f"{protein_name}_segment_{segment_index + 1}"
            sdf_path = sdf_output_dir / f"{protein_segment_id}_docking_results.sdf"

            # DockingResultCollectionのexport_to_sdfメソッドを直接呼び出す
            data["results"].export_to_sdf(sdf_path)
            # ログ出力はexport_to_sdfメソッド内で行われるため、ここでは不要

            # タンパク質構造ファイルをコピー
            try:
                copy_protein_structure(
                    protein=data["protein"], output_dir=sdf_output_dir / "proteins", name_prefix=protein_segment_id
                )
                # ログ出力はcopy_protein_structure内で行われるため、ここでは不要
            except Exception as e:
                print(f"[Protein:{protein.id}] 構造ファイルのコピー中にエラーが発生しました: {e}")
        except Exception as e:
            print(f"  SDFファイルへのエクスポート中にエラーが発生しました: {e}")

    # 統合された結果からトップヒットも表示
    combined_results = DockingResultCollection()
    for data in protein_segment_results.values():
        combined_results.extend(data["results"])

    # 統合結果のサマリーを表示
    print("\n=== 全タンパク質・セグメント統合の結果 ===")
    print(combined_results.summarize())

    # 上位ヒットの詳細情報を表示
    top_hits = combined_results.get_top(10)  # 上位10件を表示
    print(f"\n=== 全タンパク質・セグメント統合の上位 {len(top_hits)} 件 ===")

    for j, hit in enumerate(top_hits):
        # ヒットのタスクIDからタンパク質名とセグメント情報を取得
        task_id = hit.metadata.get("task_id", "unknown")
        protein_info = "不明"
        segment_info = "不明"

        # タスクIDから直接タンパク質名とセグメント情報を抽出
        if "protein_" in task_id and "segment_" in task_id:
            try:
                protein_name = task_id.split("protein_")[1].split("_segment_")[0]
                segment_str = task_id.split("_segment_")[1].split("_chunk_")[0]
                segment_index = int(segment_str)
                protein_info = protein_name
                segment_info = f"セグメント {segment_index + 1}"
            except (IndexError, ValueError):
                pass

        # 対応するグリッドボックス情報を取得
        grid_box_key = None
        if "protein_" in task_id and "segment_" in task_id:
            try:
                protein_name = task_id.split("protein_")[1].split("_segment_")[0]
                segment_index = int(task_id.split("_segment_")[1].split("_chunk_")[0])
                grid_box_key = (protein_name, segment_index)
            except (IndexError, ValueError):
                pass

        grid_box = grid_boxes.get(grid_box_key) if grid_box_key else None

        # スコア情報の表示
        scores = hit.metadata.get("scores", [])
        score_str = (
            f"全ポーズのスコア: {[float(s) for s in scores[0]]}" if scores else f"スコア: {hit.docking_score:.3f}"
        )
        print(f"   {j+1}. タンパク質: {protein_info}, {segment_info} - {score_str}")

        # グリッドボックス情報の表示
        if grid_box:
            print(f"     ドッキング範囲: {str(grid_box)}")

    # 永続化された結果の確認メッセージを表示
    print("\n=== 永続化されたドッキング結果 ===")
    print(f"リポジトリディレクトリ: {repository_dir}")
    print(f"永続化された結果は後で検索・利用できます。")

    # 合計実行時間の計算と表示
    total_execution_time = 0.0
    valid_results_count = 0
    for i, task in enumerate(docking_tasks):
        task_result = results[i]
        if hasattr(task_result, "execution_time") and task_result.execution_time is not None:
            total_execution_time += task_result.execution_time
            valid_results_count += 1

    if valid_results_count > 0:
        print("\n=== 実行時間集計 ===")
        print(f"実行時間計測対象タスク数: {valid_results_count}")
        print(f"全ドッキングタスクの合計実行時間: {total_execution_time:.2f} 秒")
        avg_time = total_execution_time / valid_results_count
        print(f"タスクあたりの平均実行時間: {avg_time:.2f} 秒")
    else:
        print("\n実行時間が計測されたタスクはありませんでした。")

    # 統合されたSDFファイルやタンパク質ごとのSDFファイルは出力せず、
    # セグメントごとの結果のみを出力します（上記のループ内で既に出力済み）
    print("\nセグメントごとのドッキング結果のみをSDFファイルにエクスポートしました。")
    try:
        # 何もしない（セグメントごとの結果は既に出力済み）
        pass
    except Exception as e:
        print(f"\nSDFファイルへのエクスポート中にエラーが発生しました: {e}")

    return results


if __name__ == "__main__":
    run_parallel_docking()
