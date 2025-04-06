"""
DaskExecutorを使って複数のドッキング計算を並列実行するサンプル。

このサンプルでは、DaskExecutorを使って複数のDockingToolABC.run_docking()
メソッド呼び出しを並列に実行する方法を示します。

Alphafoldで作成されたタンパク質構造をセグメンテーションし、
それぞれの構造に対してドッキング計算を並列実行します。
"""

import os
import tempfile
from pathlib import Path

# OpenBabelのwarningを完全に抑制するための設定
os.environ["BABEL_QUIET"] = "1"

# OpenBabelのインポート
import openbabel.pybel as pybel

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
from docking_automation.molecule import CompoundSet, Protein

# 各種ファイルのパスをハードコーディング
script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
alphafold_protein_paths = [
    script_dir / "input" / "alphafold" / "AF-A0A0A2JW93-F1-model_v4.pdb",
    script_dir / "input" / "alphafold" / "AF-A0A0H2UNG0-F1-model_v4.pdb",
]
compound_path = script_dir / "input" / "ALDR" / "actives_subset.sdf"
output_dir = script_dir / "output" / "alphafold_segmentation"


def _convert_mol2_to_sdf(mol2_path: Path) -> Path:
    """
    mol2ファイルをsdfファイルに変換する

    Args:
        mol2_path: mol2ファイルのパス

    Returns:
        変換後のsdfファイルのパス
    """
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as temp_file:
        temp_path = Path(temp_file.name)

    # mol2ファイルを読み込む
    mol = next(pybel.readfile("mol2", str(mol2_path)))

    # sdfファイルとして保存（overwrite=Trueを指定）
    mol.write("sdf", str(temp_path), overwrite=True)

    return temp_path


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
    print(f"入力タンパク質: ID={protein.id}, パス={protein.path}")

    # AlphaCutterのオプション設定
    options = {
        "loop_min": 20,
        "helix_min": 30,
        "fragment_min": 5,
        "domain_min": 50,
        "pLDDT_min": 0,
        "local_contact_range": 5,
        "domain_out": True,
        "single_out": True,
    }

    # ProteinSegmentationServiceの作成と実行
    service = ProteinSegmentationService()
    print(f"タンパク質セグメンテーションを実行中...")
    segmented_proteins = service.segment(protein, options, output_dir)

    # 結果の表示
    print(f"\nセグメンテーション完了: {len(segmented_proteins)}個のセグメントが生成されました。")
    for i, seg_protein in enumerate(segmented_proteins, 1):
        print(f"セグメント {i}: ID={seg_protein.id}, パス={seg_protein.path}")

    return segmented_proteins


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

        if not segmented_proteins:
            print(f"タンパク質 {protein_name} のセグメンテーションに失敗しました。")
            continue

        all_segmented_proteins.extend(segmented_proteins)
        protein_names.extend([protein_name] * len(segmented_proteins))

    if not all_segmented_proteins:
        print("すべてのタンパク質のセグメンテーションに失敗しました。処理を中止します。")
        return

    # 化合物セットの読み込みと分割
    chunk_size = 1  # チャンク分割の差異の、チャンクあたりの化合物数
    compound_set = CompoundSet(compound_path)
    compound_sets = compound_set.split_by_chunks(chunk_size)

    # ドッキング計算タスクの作成
    docking_tool = AutoDockVina()
    docking_tasks: list[Task] = []

    # 各セグメントに対してドッキング範囲の推定とドッキングタスクの作成
    grid_boxes = {}
    print("\n=== タンパク質・セグメントごとのドッキング範囲推定 ===")

    for i, (protein, protein_name) in enumerate(zip(all_segmented_proteins, protein_names)):
        # GridBoxの取得
        print(f"\nタンパク質: {protein_name}, セグメント {i + 1} ({protein.id}) のドッキング範囲を推定中...")
        grid_box = GridBox.from_fpocket(protein)

        # グリッドボックス情報の表示
        print(f"  中心座標: X={grid_box.center[0]:.3f}, Y={grid_box.center[1]:.3f}, Z={grid_box.center[2]:.3f}")
        print(f"  サイズ: X={grid_box.size[0]:.3f}, Y={grid_box.size[1]:.3f}, Z={grid_box.size[2]:.3f}")

        # グリッドボックスを保存
        grid_boxes[(protein_name, i)] = grid_box

        for j, split_compound_set in enumerate(compound_sets):
            # 各分割された化合物セットに対してタスクを作成
            task = Task.create(
                function=docking_tool.run_docking,
                args={
                    "protein": protein,
                    "compound_set": split_compound_set,
                    "grid_box": grid_box,
                    "additional_params": AutoDockVinaParameters(),
                },
                id=f"docking_task_protein_{protein_name}_segment_{i}_chunk_{j}",
            )
            docking_tasks.append(task)

    # タスクの登録と並列計算の実行
    executor = DaskExecutor(scheduler_type="local", n_workers=None)
    task_manager = TaskManager(executor=executor)
    for task in docking_tasks:
        task_manager.add_task(task)
    results = task_manager.execute_all()

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
        segment_top_hits = data["results"].get_top(5)  # 各セグメントの上位5件
        grid_box = grid_boxes.get((protein_name, segment_index))

        print(f"\nタンパク質: {protein_name}, セグメント {segment_index + 1} ({protein.id}) の結果:")
        if grid_box:
            print(
                f"  ドッキング範囲: 中心座標=({grid_box.center[0]:.3f}, {grid_box.center[1]:.3f}, {grid_box.center[2]:.3f}), "
                f"サイズ=({grid_box.size[0]:.3f}, {grid_box.size[1]:.3f}, {grid_box.size[2]:.3f})"
            )

        print(f"  上位 {len(segment_top_hits)} 件のドッキング結果:")
        for j, hit in enumerate(segment_top_hits):
            scores = hit.metadata["scores"]
            print(f"    全ポーズのスコア: {[float(s) for s in scores[0]]}")

    # 統合された結果からトップヒットも表示
    combined_results = DockingResultCollection()
    for data in protein_segment_results.values():
        combined_results.extend(data["results"])

    top_hits = combined_results.get_top(10)  # 上位10件を表示
    print(f"\n=== 全タンパク質・セグメント統合の上位 {len(top_hits)} 件 ===")
    for j, hit in enumerate(top_hits):
        scores = hit.metadata["scores"]
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

        print(f"   タンパク質: {protein_info}, {segment_info} - 全ポーズのスコア: {[float(s) for s in scores[0]]}")
        if grid_box:
            print(
                f"     ドッキング範囲: 中心=({grid_box.center[0]:.3f}, {grid_box.center[1]:.3f}, {grid_box.center[2]:.3f}), "
                f"サイズ=({grid_box.size[0]:.3f}, {grid_box.size[1]:.3f}, {grid_box.size[2]:.3f})"
            )

    return results


if __name__ == "__main__":
    run_parallel_docking()
