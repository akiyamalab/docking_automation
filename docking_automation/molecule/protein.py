from pathlib import Path
from typing import Optional


class Protein:
    """
    タンパク質を表すクラス。
    現時点では .pdb 形式をサポートしている。
    """
    def __init__(self, path: str|Path, id: Optional[str] = None):
        """
        Proteinオブジェクトを初期化する。
        
        Args:
            path: タンパク質のファイルパス
            id: タンパク質のID（指定しない場合はファイル名がIDとなる）
        """
        self.path: Path = Path(path)
        self.id: str = id if id is not None else self.path.stem