import gzip
import hashlib
import io
import os
from pathlib import Path
from typing import Iterator, List, Tuple, Union


def expandpath(path: str | Path) -> Path:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = Path(path)
    path = path.expanduser()
    return Path(os.path.expandvars(path))


def safe_open(filepath: str | Path, mode: str = "rt") -> Union[gzip.GzipFile, io.TextIOWrapper]:
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
        fileobj = open(filepath, mode)  # type: ignore
    return fileobj


def read_compounds_from_sdf(filepath: str | Path) -> Iterator[Tuple[int, List[str | bytes]]]:
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


def calculate_file_content_hash(filepath: str | Path) -> str:
    """
    ファイルの内容に基づいたハッシュ値を計算する。

    SHA-256アルゴリズムを使用して、ファイルの内容からハッシュ値を生成する。

    Args:
        filepath: ハッシュ値を計算するファイルのパス

    Returns:
        str: ファイル内容のSHA-256ハッシュ値（16進数文字列）

    Raises:
        ValueError: ファイルが存在しない場合
    """
    filepath = expandpath(filepath)
    if not Path(filepath).exists():
        raise ValueError(f"指定されたパス '{filepath}' が存在しません")

    sha256_hash = hashlib.sha256()

    # ファイルを適切なモードで開く（バイナリモード）
    with open(filepath, "rb") as f:
        # 一度に読み込むバイト数を指定（大きなファイルでもメモリ効率が良い）
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)

    return sha256_hash.hexdigest()
