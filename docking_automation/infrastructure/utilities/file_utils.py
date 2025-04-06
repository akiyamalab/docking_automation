import gzip
import os
from pathlib import Path


def expandpath(path: str | Path) -> Path:
    """
    Expand ~ and $HOME and other environment variables
    """
    path = Path(path)
    path = path.expanduser()
    return Path(os.path.expandvars(path))


def open(filepath: str | Path):
    filepath = expandpath(filepath)
    if filepath.suffix == ".gz":
        fileobj = gzip.open(filepath, "rt")
    else:
        fileobj = open(filepath)
    return fileobj
