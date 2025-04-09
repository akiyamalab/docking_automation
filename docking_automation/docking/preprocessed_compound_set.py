from pathlib import Path
from typing import Any, Dict, List, Optional


# 値オブジェクト
class PreprocessedCompoundSet:
    """
    前処理済みの化合物群を表すオブジェクト

    複数の化合物ファイルパスを保持できるようにしている。
    ファイルベース、メモリベースどちらでも対応できるようにしている。
    実際にどのフィールドを使うかはドッキングツールによって依存する。
    """

    def __init__(self, file_paths: Optional[List[Path]] = None, data: Optional[object] = None):
        """
        前処理済みの化合物群を表すオブジェクトを初期化する。

        Args:
            file_paths: 前処理済み化合物ファイルのパスのリスト
            data: メモリ上のデータ（オプション）
        """
        self.file_paths = file_paths or []  # 複数の化合物のファイルパス
        self.data = data

    def get_properties(self) -> Dict[str, Any]:
        """
        前処理済み化合物セットのプロパティを取得する。

        Returns:
            プロパティの辞書
        """
        # 基本的なプロパティを返す
        properties = {
            "file_paths": [str(path) for path in self.file_paths],
            "file_count": len(self.file_paths),
        }

        # データがある場合は、その情報も含める
        if self.data is not None:
            properties["has_data"] = True

        return properties
