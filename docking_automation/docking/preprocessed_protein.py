from functools import cached_property
from pathlib import Path
from typing import Optional

from docking_automation.infrastructure.utilities.file_utils import (
    calculate_file_content_hash,
)


# 値オブジェクト
class PreprocessedProtein:
    """
    前処理済みのタンパク質を表すオブジェクト

    ファイルベース、メモリベースどちらでも対応できるようにしている。
    実際にどのフィールドを使うかはドッキングツールによって依存する。
    """

    def __init__(self, file_path: Path = None, data: Optional[object] = None):
        """
        前処理済みのタンパク質を表すオブジェクトを初期化する。

        Args:
            file_path: 前処理済みタンパク質ファイルのパス
            data: メモリ上のデータ（オプション）
        """
        self.file_path = file_path
        self.data = data

    @cached_property
    def content_hash(self) -> str:
        """
        ファイルの内容に基づいたハッシュ値を取得する。

        初回アクセス時に計算し、その後はキャッシュした値を返す。

        Returns:
            ファイル内容のSHA-256ハッシュ値（16進数文字列）
        """
        return calculate_file_content_hash(self.file_path)
