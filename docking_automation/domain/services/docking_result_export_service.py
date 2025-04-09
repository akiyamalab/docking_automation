"""
ドッキング結果のエクスポートサービスモジュール。

このモジュールは、ドッキング結果をSDFファイルにエクスポートするためのドメインサービスを提供します。
"""

from pathlib import Path
from typing import Optional

from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepositoryABC,
)


def export_docking_results_for_protein(
    repository: DockingResultRepositoryABC, protein_id: str, output_path: Optional[Path] = None
) -> Path:
    """
    指定されたタンパク質に対するドッキング結果をSDFファイルにエクスポートする。

    Args:
        repository: ドッキング結果リポジトリ
        protein_id: タンパク質ID
        output_path: 出力ファイルのパス（指定しない場合はデフォルトのパスを使用）

    Returns:
        エクスポートされたSDFファイルのパス

    Raises:
        ValueError: タンパク質IDに対するドッキング結果が見つからない場合
    """
    # リポジトリからタンパク質IDに基づいてドッキング結果を取得
    results = repository.find_by_protein(protein_id)

    # 結果が空の場合はエラー
    if len(results) == 0:
        raise ValueError(f"タンパク質ID '{protein_id}' に対するドッキング結果が見つかりません")

    # 出力パスが指定されていない場合はデフォルトのパスを使用
    if output_path is None:
        output_path = Path(f"{protein_id}_docking_results.sdf")

    # 出力ディレクトリが存在しない場合は作成
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # DockingResultCollectionのexport_to_sdfメソッドを呼び出す
    results.export_to_sdf(output_path)

    return output_path
