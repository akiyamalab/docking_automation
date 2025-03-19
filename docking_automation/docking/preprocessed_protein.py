from pathlib import Path
from typing import Optional, List

# 値オブジェクト
class PreprocessedProtein:
    """
    前処理済みのタンパク質を表すオブジェクト
    
    ファイルベース、メモリベースどちらでも対応できるようにしている。
    実際にどのフィールドを使うかはドッキングツールによって依存する。
    """
    
    def __init__(self, file_paths: Optional[List[Path]] = None, data: Optional[object] = None):
        """
        前処理済みのタンパク質を表すオブジェクトを初期化する。
        
        Args:
            file_paths: 前処理済みタンパク質ファイルのパスのリスト
            data: メモリ上のデータ（オプション）
        """
        self.file_paths = file_paths or []  # 複数のタンパク質ファイルのパス
        self.data = data