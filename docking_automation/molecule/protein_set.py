from __future__ import annotations

import uuid
from pathlib import Path
from typing import Callable, Dict, Iterable, Iterator, List, Optional

from docking_automation.molecule.protein import Protein


def _count_ca_atoms(pdb_path: Path) -> int:
    """PDB ファイルから CA 原子数をカウント（残基数の近似）。"""
    count = 0
    with open(pdb_path) as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                count += 1
    return count


class ProteinSet:
    """複数 Protein の集約ルート。

    不変条件:
        - protein_id はセット内で一意
        - 各 Protein の path は実在ファイル
    """

    def __init__(self, proteins: Iterable[Protein], id: Optional[str] = None) -> None:
        protein_list = list(proteins)
        seen: Dict[str, Protein] = {}
        for p in protein_list:
            if p.id in seen:
                raise ValueError(f"protein_id '{p.id}' が重複しています")
            seen[p.id] = p
        self._proteins: Dict[str, Protein] = seen
        self._order: List[str] = [p.id for p in protein_list]
        self._id: str = id if id is not None else str(uuid.uuid4())
        self._content_hashes_cache: Optional[Dict[str, str]] = None

    def __iter__(self) -> Iterator[Protein]:
        return (self._proteins[pid] for pid in self._order)

    def __len__(self) -> int:
        return len(self._order)

    def __contains__(self, protein_id: object) -> bool:
        return protein_id in self._proteins

    def __getitem__(self, protein_id: str) -> Protein:
        if protein_id not in self._proteins:
            raise KeyError(f"protein_id '{protein_id}' が見つかりません")
        return self._proteins[protein_id]

    @property
    def id(self) -> str:
        return self._id

    @property
    def protein_ids(self) -> List[str]:
        return list(self._order)

    def content_hashes(self) -> Dict[str, str]:
        if self._content_hashes_cache is None:
            self._content_hashes_cache = {
                pid: self._proteins[pid].content_hash for pid in self._order
            }
        return dict(self._content_hashes_cache)

    def subset(self, protein_ids: Iterable[str]) -> "ProteinSet":
        ids = list(protein_ids)
        for pid in ids:
            if pid not in self._proteins:
                raise KeyError(f"protein_id '{pid}' が見つかりません")
        return ProteinSet([self._proteins[pid] for pid in ids])

    def filter(self, predicate: Callable[[Protein], bool]) -> "ProteinSet":
        return ProteinSet([p for p in self if predicate(p)])

    @classmethod
    def from_directory(
        cls,
        path: str | Path,
        pattern: str = "*.pdb",
        recursive: bool = False,
        id: Optional[str] = None,
    ) -> "ProteinSet":
        root = Path(path)
        files = sorted(root.rglob(pattern) if recursive else root.glob(pattern))
        proteins = [Protein(f) for f in files]
        return cls(proteins, id=id)

    @classmethod
    def from_afdb_mouse(
        cls,
        root: str | Path = "data/afdb/pdb",
        limit: Optional[int] = None,
        max_residues: Optional[int] = 2000,
        id: str = "afdb_mouse",
    ) -> "ProteinSet":
        root_path = Path(root)
        skipped_path = root_path.parent / "skipped_large.txt"

        files = sorted(root_path.glob("*.pdb"))
        proteins: List[Protein] = []
        skipped: List[str] = []

        for pdb_file in files:
            if limit is not None and len(proteins) >= limit:
                break
            if max_residues is not None:
                ca_count = _count_ca_atoms(pdb_file)
                if ca_count > max_residues:
                    skipped.append(pdb_file.name)
                    continue
            proteins.append(Protein(pdb_file))

        if skipped:
            skipped_path.parent.mkdir(parents=True, exist_ok=True)
            with open(skipped_path, "a") as f:
                for name in skipped:
                    f.write(name + "\n")

        return cls(proteins, id=id)
