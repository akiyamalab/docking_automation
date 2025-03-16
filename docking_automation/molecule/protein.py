from pathlib import Path


class Protein:
    """
    タンパク質を表すクラス。
    現時点では .pdb 形式をサポートしている。
    """
    def __init__(self, path: str|Path):
        self.path: Path = Path(path)