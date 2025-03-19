from pathlib import Path
from typing import Optional, List

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