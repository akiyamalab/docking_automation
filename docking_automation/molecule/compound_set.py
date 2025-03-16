from pathlib import Path


class CompoundSet:
    """
    化合物を表すクラス。
    現時点では .sdf 形式をサポートしている。
    """
    def __init__(self, path: str|Path):
        self.path: Path = Path(path)