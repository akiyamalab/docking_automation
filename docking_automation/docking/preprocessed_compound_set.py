import hashlib
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, List, Optional

from docking_automation.infrastructure.utilities.file_utils import (
    calculate_file_content_hash,
)


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

    @cached_property
    def content_hash(self) -> str:
        """
        ファイルの内容に基づいたハッシュ値を取得する。

        初回アクセス時に計算し、その後はキャッシュした値を返す。

        Returns:
            ファイル内容のSHA-256ハッシュ値（16進数文字列）
        """
        hashes = [calculate_file_content_hash(file_path) for file_path in self.file_paths]
        # 複数のファイルのハッシュを結合して最終的なハッシュを計算
        combined_hash = "".join(hashes)
        sha256_hash = hashlib.sha256()
        sha256_hash.update(combined_hash.encode("utf-8"))
        combined_hash = sha256_hash.hexdigest()
        return combined_hash
