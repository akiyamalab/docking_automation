import gzip
from pathlib import Path
from typing import Optional, TextIO, Union, BinaryIO
from rdkit import Chem


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
        # ファイル形式を判断
        file_format = self.path.suffix.lower()
        
        # gzipで圧縮されている場合
        if file_format == '.gz':
            # 実際のファイル形式を取得
            actual_format = self.path.stem.split('.')[-1].lower()
            
            if actual_format == 'sdf':
                return self._count_compounds_in_sdf()
            else:
                raise ValueError(f"サポートされていないファイル形式です: {actual_format}")
        # 非圧縮ファイルの場合
        elif file_format == '.sdf':
            return self._count_compounds_in_sdf()
        else:
            raise ValueError(f"サポートされていないファイル形式です: {file_format}")
    
    def _count_compounds_in_sdf(self) -> int:
        """
        SDFファイル内の化合物数をカウントする。
        RDKitを使用して化合物数をカウントする。

        Returns:
            化合物の数
        """
        try:
            if str(self.path).endswith('.gz'):
                with gzip.open(self.path, 'rt') as f:
                    count = 0
                    for line in f:
                        if "$$$$" in line:
                            count += 1
                    return count
            else:
                with open(self.path, 'rt') as f:
                    count = 0
                    for line in f:
                        if "$$$$" in line:
                            count += 1
                    return count
        except Exception as e:
            raise ValueError(f"化合物数のカウント中にエラーが発生しました: {e}")