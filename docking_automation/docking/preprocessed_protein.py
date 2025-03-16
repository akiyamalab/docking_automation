from pathlib import Path
from typing import Optional


class PreprocessedProtein:
    """
    前処理済みのタンパク質を表すオブジェクト
    ファイルベース、メモリベースどちらでも対応できるようにしている。実際にどのフィールドを使うかはドッキングツールによって依存する。
    """
    
    def __init__(self, file_path: Optional[Path] = None, data: Optional[object] = None):
        """
        前処理済みのタンパク質を表すオブジェクト
        ファイルベース、メモリベースどちらでも対応できるようにしている。実際にどのフィールドを使うかはドッキングツールによって依存する。
        """
        self.file_path = file_path
        self.data = data