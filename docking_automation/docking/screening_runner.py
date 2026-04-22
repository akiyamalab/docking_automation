from __future__ import annotations

import gzip as gz
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterator, List, TYPE_CHECKING

if TYPE_CHECKING:
    from docking_automation.molecule.protein_set import ProteinSet
    from docking_automation.molecule.compound_set import CompoundSet
    from docking_automation.docking.grid_box_cache import GridBoxCache
    from docking_automation.molecule.protein import Protein
    from docking_automation.docking.grid_box import GridBox


@dataclass(frozen=True)
class ScreeningResult:
    total_pairs: int
    new_pairs: int
    reused_pairs: int
    failed_pairs: int
    elapsed_sec: float
    hdf5_path: Path
    log_path: Path


class ScreeningRunner:
    """N×M ドッキング司令塔。再開可能・冪等。Wave 1: 逐次実行スタブ。"""

    def __init__(
        self,
        protein_set: "ProteinSet",
        compound_set: "CompoundSet",
        grid_box_cache: "GridBoxCache",
        hdf5_path: Path,
        dask_n_workers: int = 4,
        exhaustiveness: int = 1,
        log_path: Path = Path("logs/screening_run.jsonl"),
        save_poses: bool = True,
        top_n_poses: int = 1,
        compression: str = "gzip",
        grid_box_missing_policy: str = "skip",
    ) -> None:
        self.protein_set = protein_set
        self.compound_set = compound_set
        self.grid_box_cache = grid_box_cache
        self.hdf5_path = Path(hdf5_path)
        self.dask_n_workers = dask_n_workers
        self.exhaustiveness = exhaustiveness
        self.log_path = Path(log_path)
        self.save_poses = save_poses
        self.top_n_poses = top_n_poses
        self.compression = compression
        self.grid_box_missing_policy = grid_box_missing_policy

    def run(self, resume: bool = True, _repo: Any = None) -> ScreeningResult:
        """実行メインループ。Wave 1: 逐次実行。_repo はテスト用依存注入。"""
        t0 = time.monotonic()
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

        if not resume and self.hdf5_path.exists():
            self.hdf5_path.unlink()

        if _repo is None:
            from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
                HDF5DockingResultRepository,
            )
            _repo = HDF5DockingResultRepository(self.hdf5_path, mode="append")
        repo = _repo
        protein_hashes = self.protein_set.content_hashes()

        all_pairs = list(self._enumerate_pairs())
        total = len(all_pairs)
        unprocessed = self._filter_unprocessed(repo)
        reused = total - len(unprocessed)
        new_pairs = 0
        failed = 0

        with open(self.log_path, "a") as log_fp:
            for protein_id, compound_index in unprocessed:
                protein = self.protein_set[protein_id]
                grid_box = self.grid_box_cache.get(protein)
                compound_hash = self.compound_set.get_compound_hash(compound_index)
                ts = datetime.now(timezone.utc).isoformat()

                if grid_box is None:
                    if self.grid_box_missing_policy == "error":
                        raise ValueError(f"GridBox missing for protein '{protein_id}'")
                    log_fp.write(
                        json.dumps({
                            "ts": ts,
                            "protein_id": protein_id,
                            "compound_idx": compound_index,
                            "compound_hash": compound_hash,
                            "error": "grid_box_missing",
                        }) + "\n"
                    )
                    continue

                try:
                    score, pose_blob = self._dock_one_pair(protein, compound_index, grid_box)
                    p_hash = protein_hashes[protein_id]
                    self._save_to_hdf5(p_hash, protein_id, compound_index, compound_hash, score, pose_blob)
                    new_pairs += 1
                    log_fp.write(
                        json.dumps({
                            "ts": ts,
                            "protein_id": protein_id,
                            "compound_idx": compound_index,
                            "compound_hash": compound_hash,
                            "score": float(score),
                            "status": "new",
                        }) + "\n"
                    )
                except Exception as e:
                    failed += 1
                    log_fp.write(
                        json.dumps({
                            "ts": ts,
                            "protein_id": protein_id,
                            "compound_idx": compound_index,
                            "compound_hash": compound_hash,
                            "error": str(e),
                            "status": "failed",
                        }) + "\n"
                    )

        elapsed = time.monotonic() - t0
        return ScreeningResult(
            total_pairs=total,
            new_pairs=new_pairs,
            reused_pairs=reused,
            failed_pairs=failed,
            elapsed_sec=elapsed,
            hdf5_path=self.hdf5_path,
            log_path=self.log_path,
        )

    def _enumerate_pairs(self) -> Iterator[tuple[str, int]]:
        """(protein_id, compound_index) の全組合せを yield"""
        for protein in self.protein_set:
            for i in range(self.compound_set.get_compound_count()):
                yield (protein.id, i)

    def _filter_unprocessed(self, repo: Any) -> List[tuple[str, int]]:
        """HDF5 既存キーでフィルタし、未実行ペアのみ返す。"""
        protein_hashes = self.protein_set.content_hashes()
        result = []
        for protein_id, compound_index in self._enumerate_pairs():
            p_hash = protein_hashes[protein_id]
            c_hash = self.compound_set.get_compound_hash(compound_index)
            if not repo._exists(p_hash, c_hash):
                result.append((protein_id, compound_index))
        return result

    def _dock_one_pair(
        self,
        protein: "Protein",
        compound_index: int,
        grid_box: "GridBox",
    ) -> tuple[float, bytes]:
        """Wave 1: 逐次 Vina ドッキング。Wave 2 で Dask 対応に差し替える。"""
        import tempfile

        from vina import Vina

        from docking_automation.converters.molecule_converter import MoleculeConverter
        from docking_automation.infrastructure.utilities.file_utils import read_compounds_from_sdf

        converter = MoleculeConverter()
        temp_dir = Path(tempfile.mkdtemp())

        pdbqt_path = temp_dir / f"{protein.id}.pdbqt"
        converter.protein_to_pdbqt(protein, pdbqt_path)

        compound_sdf = temp_dir / f"compound_{compound_index}.sdf"
        for i, (_, lines) in enumerate(read_compounds_from_sdf(self.compound_set.path)):
            if i == compound_index:
                compound_sdf.write_text("".join(str(l) for l in lines))
                break

        compound_pdbqt = temp_dir / f"compound_{compound_index}.pdbqt"
        converter.sdf_to_pdbqt(compound_sdf, compound_pdbqt)

        center = grid_box.center
        size = grid_box.size

        v = Vina(cpu=1, seed=1, verbosity=0)
        v.set_receptor(str(pdbqt_path))
        v.set_ligand_from_file(str(compound_pdbqt))
        v.compute_vina_maps(
            center=[center[0], center[1], center[2]],
            box_size=[size[0], size[1], size[2]],
        )
        v.dock(exhaustiveness=self.exhaustiveness, n_poses=self.top_n_poses, min_rmsd=1.0)

        output_pdbqt = temp_dir / f"output_{compound_index}.pdbqt"
        output_sdf = temp_dir / f"output_{compound_index}.sdf"
        v.write_poses(str(output_pdbqt), n_poses=self.top_n_poses, overwrite=True)
        converter.pdbqt_to_sdf(output_pdbqt, output_sdf)

        scores = v.energies()
        score = float(scores[0, 0])
        pose_bytes = output_sdf.read_bytes()
        pose_blob = gz.compress(pose_bytes, compresslevel=9)
        return score, pose_blob

    def _save_to_hdf5(
        self,
        protein_hash: str,
        protein_id: str,
        compound_index: int,
        compound_hash: str,
        score: float,
        pose_blob: bytes,
    ) -> None:
        """Phase 2 スキーマで HDF5 に保存する。"""
        import h5py

        self.hdf5_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(self.hdf5_path, "a", libver="latest") as f:
            group_path = f"/results/{protein_hash}/{compound_hash}"
            if group_path in f:
                return
            g = f.require_group(group_path)
            g.attrs["protein_id"] = protein_id
            g.attrs["compound_index"] = compound_index
            g.create_dataset("docking_score", data=float(score), dtype="f4")
            g.create_dataset("pose_blob", data=pose_blob)
            g.create_dataset(
                "computed_at",
                data=datetime.now(timezone.utc).isoformat(),
                dtype=h5py.string_dtype(),
            )
            g.create_dataset("source", data="vina", dtype=h5py.string_dtype())
            g.create_dataset("top_n", data=self.top_n_poses, dtype="i4")
