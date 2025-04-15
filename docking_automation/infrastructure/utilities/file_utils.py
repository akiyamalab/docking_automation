import gzip
import os
from pathlib import Path
from typing import Iterator, List, Tuple


def expandpath(path: str | Path) -> Path:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = Path(path)
    path = path.expanduser()
    return Path(os.path.expandvars(path))


def safe_open(filepath: str | Path, mode: str = "rt"):
    """
    ファイルを安全に開くためのユーティリティ関数。
    gzipで圧縮されたファイルの場合は自動的に展開して開きます。

    Args:
        filepath: 開くファイルのパス
        mode: ファイルを開くモード（デフォルトは"rt"）

    Returns:
        開いたファイルオブジェクト
    """
    filepath = expandpath(filepath)
    if str(filepath).endswith(".gz"):
        fileobj = gzip.open(filepath, mode)
    else:
        fileobj = __builtins__["open"](filepath, mode)
    return fileobj


def read_compounds_from_sdf(filepath: str | Path) -> Iterator[Tuple[int, List[str]]]:
    """
    SDFファイルから化合物データを読み込み、各化合物データをyieldするジェネレータ。
    各化合物は"$$$$"で区切られる。

    Args:
        filepath: SDFファイルのパス

    Yields:
        Tuple[int, List[str]]: (化合物のインデックス, 化合物の行のリスト)

    Raises:
        ValueError: ファイルの読み込み中にエラーが発生した場合
    """
    try:
        with safe_open(filepath) as f:
            compound_index = 0
            current_compound = []

            for line in f:
                current_compound.append(line)
                if "$$$$" in line:
                    # 化合物の終端に達した場合
                    yield compound_index, current_compound
                    current_compound = []
                    compound_index += 1
    except Exception as e:
        raise ValueError(f"SDFファイルの読み込み中にエラーが発生しました: {e}")


def count_compounds_in_sdf(filepath: str | Path) -> int:
    """
    SDFファイル内の化合物数をカウントする。

    Args:
        filepath: SDFファイルのパス

    Returns:
        int: 化合物の数

    Raises:
        ValueError: 化合物数のカウント中にエラーが発生した場合
    """
    try:
        count = sum(1 for _ in read_compounds_from_sdf(filepath))
        return count
    except Exception as e:
        raise ValueError(f"化合物数のカウント中にエラーが発生しました: {e}")
