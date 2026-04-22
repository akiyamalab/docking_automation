from __future__ import annotations

import json
import os
import tempfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Dict, Iterator, List, Optional

from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.protein import Protein

if TYPE_CHECKING:
    from docking_automation.molecule.protein_set import ProteinSet


@dataclass(frozen=True)
class GridBoxCacheEntry:
    protein_id: str
    protein_content_hash: str
    grid_box: GridBox
    source: str
    computed_at: str

    def to_dict(self) -> dict:
        return {
            "protein_content_hash": self.protein_content_hash,
            "center": self.grid_box.center.tolist(),
            "size": self.grid_box.size.tolist(),
            "source": self.source,
            "computed_at": self.computed_at,
        }

    @classmethod
    def from_dict(cls, protein_id: str, d: dict) -> "GridBoxCacheEntry":
        grid_box = GridBox(center=tuple(d["center"]), size=tuple(d["size"]))
        return cls(
            protein_id=protein_id,
            protein_content_hash=d["protein_content_hash"],
            grid_box=grid_box,
            source=d["source"],
            computed_at=d["computed_at"],
        )


class GridBoxCache:
    def __init__(
        self,
        path: str | Path,
        predictor: str = "fpocket",
        predictor_version: str = "4.0",
        pocket_rank: int = 1,
    ) -> None:
        self._path = Path(path)
        self._predictor = predictor
        self._predictor_version = predictor_version
        self._pocket_rank = pocket_rank
        self._created_at = datetime.now(timezone.utc).isoformat()
        self._entries: Dict[str, GridBoxCacheEntry] = {}

    def get(self, protein: Protein) -> Optional[GridBox]:
        entry = self._entries.get(protein.id)
        if entry is None:
            return None
        if entry.protein_content_hash != protein.content_hash:
            return None
        return entry.grid_box

    def has(self, protein: Protein) -> bool:
        entry = self._entries.get(protein.id)
        if entry is None:
            return False
        return entry.protein_content_hash == protein.content_hash

    def put(self, protein: Protein, grid_box: GridBox, source: str = "fpocket") -> None:
        computed_at = datetime.now(timezone.utc).isoformat()
        self._entries[protein.id] = GridBoxCacheEntry(
            protein_id=protein.id,
            protein_content_hash=protein.content_hash,
            grid_box=grid_box,
            source=source,
            computed_at=computed_at,
        )

    def invalidate(self, protein_id: str) -> None:
        self._entries.pop(protein_id, None)

    def __len__(self) -> int:
        return len(self._entries)

    def __contains__(self, protein_id: str) -> bool:
        return protein_id in self._entries

    def missing_ids(self, protein_set: "ProteinSet") -> List[str]:
        result = []
        for protein in protein_set:
            entry = self._entries.get(protein.id)
            if entry is None or entry.protein_content_hash != protein.content_hash:
                result.append(protein.id)
        return result

    def entries(self) -> Iterator[GridBoxCacheEntry]:
        return iter(self._entries.values())

    def _to_dict(self) -> dict:
        return {
            "version": 1,
            "created_at": self._created_at,
            "predictor": self._predictor,
            "predictor_version": self._predictor_version,
            "pocket_rank": self._pocket_rank,
            "entries": {pid: entry.to_dict() for pid, entry in self._entries.items()},
        }

    def save(self, path: Optional[str | Path] = None, atomic: bool = True) -> None:
        target = Path(path or self._path)
        target.parent.mkdir(parents=True, exist_ok=True)
        data = self._to_dict()
        if atomic:
            with tempfile.NamedTemporaryFile(
                mode="w", dir=target.parent, suffix=".tmp", delete=False
            ) as f:
                json.dump(data, f, indent=2)
                f.flush()
                os.fsync(f.fileno())
                tmp_path = f.name
            os.replace(tmp_path, target)
        else:
            with open(target, "w") as f:
                json.dump(data, f, indent=2)

    @classmethod
    def from_file(cls, path: str | Path) -> "GridBoxCache":
        p = Path(path)
        cache = cls(path=p)
        if not p.exists():
            return cache
        with open(p) as f:
            data = json.load(f)
        cache._predictor = data.get("predictor", "fpocket")
        cache._predictor_version = data.get("predictor_version", "4.0")
        cache._pocket_rank = data.get("pocket_rank", 1)
        cache._created_at = data.get("created_at", cache._created_at)
        for protein_id, entry_dict in data.get("entries", {}).items():
            cache._entries[protein_id] = GridBoxCacheEntry.from_dict(protein_id, entry_dict)
        return cache

    @classmethod
    def build_from_protein_set(
        cls,
        protein_set: "ProteinSet",
        predictor,
        cache_path: str | Path,
        save_every: int = 100,
        on_error: str = "skip",
    ) -> "GridBoxCache":
        raise NotImplementedError("build_from_protein_set は Phase 2 で実装予定")
