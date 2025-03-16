from pathlib import Path
from typing import Optional


class CompoundSet:
    """
    化合物セットを表すクラス。
    現時点では .sdf 形式をサポートしている。
    1つのファイルに複数の化合物が含まれる。
    """
    def __init__(self, path: str|Path, id: Optional[str] = None):
        """
        CompoundSetオブジェクトを初期化する。
        
        Args:
            path: 化合物セットのファイルパス
            id: 化合物セットのID（指定しない場合はファイル名がIDとなる）
        """
        self.path: Path = Path(path)
        self.id: str = id if id is not None else self.path.stem
    
    def get_compound_count(self) -> int:
        """
        化合物セットに含まれる化合物の数を取得する。
        
        Returns:
            化合物の数
        """
        raise NotImplementedError()