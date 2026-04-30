from __future__ import annotations

import gzip as gz
import json
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, TYPE_CHECKING

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


def _make_docking_tool(backend: str):
    """backendに応じて ScreeningTool 派生クラスを生成するファクトリ。

    Returns:
        `ScreeningTool` instance (AutoDockVina / UniDockDocking / UniDock2Docking)。
    """
    if backend == "vina":
        from docking_automation.docking.autodockvina_docking import AutoDockVina
        return AutoDockVina()
    elif backend == "unidock":
        from docking_automation.docking.unidock_docking import UniDockDocking
        return UniDockDocking()
    elif backend == "unidock2":
        from docking_automation.docking.unidock2_docking import UniDock2Docking
        return UniDock2Docking()
    else:
        raise ValueError(f"Unknown backend: {backend}")


def _prepare_cache_for_dask(
    protein_path: str,
    protein_content_hash: str,
    grid_center: List[float],
    grid_size: List[float],
    cache_dir: str,
    force: bool = False,
) -> dict:
    """Dask worker で 1 受容体の UniDock2 receptor cache を生成する。

    Returns:
        {"protein_content_hash": str, "cache_path": str|None, "error": str|None}
    """
    import os
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    try:
        from docking_automation.docking.grid_box import GridBox
        from docking_automation.docking.unidock2_docking import UniDock2Docking
        from docking_automation.molecule.protein import Protein

        protein = Protein(Path(protein_path))
        grid_box = GridBox(
            center=(float(grid_center[0]), float(grid_center[1]), float(grid_center[2])),
            size=(float(grid_size[0]), float(grid_size[1]), float(grid_size[2])),
        )
        cache_path = Path(cache_dir) / f"{protein_content_hash}.json"
        tool = UniDock2Docking(cache_dir=Path(cache_dir))
        tool.prepare_receptor_cache(protein, grid_box, cache_path, force=force)
        return {
            "protein_content_hash": protein_content_hash,
            "cache_path": str(cache_path),
            "error": None,
        }
    except Exception as e:
        return {
            "protein_content_hash": protein_content_hash,
            "cache_path": None,
            "error": f"{type(e).__name__}: {e}",
        }


def dock_one_protein_via_tool(
    protein_path: str,
    protein_id: str,
    protein_content_hash: str,
    compound_sdf_path: str,
    compound_indices: List[int],
    compound_hashes: Dict[int, str],
    grid_center: List[float],
    grid_size: List[float],
    backend: str = "vina",
    **kwargs: Any,
) -> List[dict]:
    """ScreeningTool 経由で 1 タンパク質ドッキングを Dask worker 上で実行する。

    ScreeningTool 経由で全 backend (vina / unidock / unidock2) に対応。
    ScreeningRunner のデフォルト dock 関数。

    Args:
        backend: "vina" / "unidock" / "unidock2" のいずれか。
        **kwargs: cache_dir (str|None), exhaustiveness, rescue_mode 等。
    """
    import gzip as gz
    import tempfile
    import time
    from pathlib import Path

    if not compound_indices:
        return []

    from docking_automation.docking.grid_box import GridBox
    from docking_automation.infrastructure.utilities.file_utils import read_compounds_from_sdf
    from docking_automation.molecule.protein import Protein

    converter = None
    if backend != 'unidock2':
        from docking_automation.converters.molecule_converter import MoleculeConverter
        converter = MoleculeConverter()
    tool = _make_docking_tool(backend)
    grid_box = GridBox(
        center=(float(grid_center[0]), float(grid_center[1]), float(grid_center[2])),
        size=(float(grid_size[0]), float(grid_size[1]), float(grid_size[2])),
    )
    protein_obj = Protein(Path(protein_path), id=protein_id)

    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        # ligand を backend 固有形式に変換 (Vina/v1: PDBQT, v2: SDF)
        ligand_paths: List[Path] = []
        compound_to_hash: Dict[Path, str] = {}
        compound_to_idx: Dict[Path, int] = {}
        idx_set = set(compound_indices)
        for i, (_, lines) in enumerate(read_compounds_from_sdf(Path(compound_sdf_path))):
            if i not in idx_set:
                continue
            sdf_path = tmp_dir / f'compound_{i}.sdf'
            sdf_path.write_text(''.join(str(l) for l in lines))
            if backend == 'unidock2':
                ligand_paths.append(sdf_path)
                target = sdf_path
            else:
                assert converter is not None  # backend != 'unidock2' なので必ず初期化済
                pdbqt_path = tmp_dir / f'compound_{i}.pdbqt'
                try:
                    converter.sdf_to_pdbqt(sdf_path, pdbqt_path)
                    ligand_paths.append(pdbqt_path)
                    target = pdbqt_path
                except Exception:
                    continue
            compound_to_hash[target] = compound_hashes.get(i, sdf_path.stem)
            compound_to_idx[target] = i

        # cache_dir が指定されていれば pre-generated cache を使う (unidock2 のみ)
        ext_cache_dir = kwargs.get('cache_dir')
        if ext_cache_dir and backend == 'unidock2':
            cache_prefix = Path(ext_cache_dir) / f'{protein_content_hash}.json'
        else:
            cache_ext = '.json' if backend == 'unidock2' else ''
            cache_prefix = tmp_dir / f'{protein_content_hash}_cache{cache_ext}'

        t0 = time.monotonic()
        if not cache_prefix.exists():
            tool.prepare_receptor_cache(protein_obj, grid_box, cache_prefix)
        dock_results = tool.dock_with_cache(
            cache=cache_prefix,
            ligand_paths=ligand_paths,
            grid_box=grid_box,
            protein_content_hash=protein_content_hash,
            compound_content_hashes=[compound_to_hash[p] for p in ligand_paths],
        )
        elapsed_per = round((time.monotonic() - t0) / max(len(ligand_paths), 1), 3)

    # 結果を ScreeningRunner._collect_and_save が期待する dict 形式に変換
    results_by_idx: Dict[int, dict] = {}
    for r in dock_results:
        lig_path = ligand_paths[r.compound_index] if r.compound_index < len(ligand_paths) else None
        orig_idx = compound_to_idx.get(lig_path) if lig_path else None
        if orig_idx is None:
            continue
        pose_blob = None
        if r.result_path and r.result_path.exists():
            pose_blob = gz.compress(r.result_path.read_bytes(), compresslevel=9)
        results_by_idx[orig_idx] = {
            'protein_id': protein_id,
            'compound_index': orig_idx,
            'protein_content_hash': protein_content_hash,
            'compound_content_hash': compound_hashes.get(orig_idx),
            'score': r.docking_score,
            'pose_blob': pose_blob,
            'elapsed_sec': elapsed_per,
            'error': None,
        }

    return [
        results_by_idx.get(idx, {
            'protein_id': protein_id,
            'compound_index': idx,
            'protein_content_hash': protein_content_hash,
            'compound_content_hash': compound_hashes.get(idx),
            'score': None,
            'pose_blob': None,
            'elapsed_sec': elapsed_per if 'elapsed_per' in dir() else 0.0,
            'error': 'not_in_dock_results',
        })
        for idx in compound_indices
    ]


class ScreeningRunner:
    """N×M ドッキング司令塔。再開可能・冪等。Wave 2: Dask LocalCluster並列実行。"""

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
        backend: str = "vina",
        search_mode: str = "balance",
        schema_version: str = "v2",
        extra_padding: float = 5.0,
        rescue_mode: bool = False,
        cache_dir: Optional[Path] = None,
        force_cache: bool = False,
        cache_n_workers: Optional[int] = None,
        _dock_fn: Optional[Callable] = None,
        _cluster_kwargs: Optional[dict] = None,
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
        self.backend = backend
        self.search_mode = search_mode
        self.schema_version = schema_version
        self.extra_padding = extra_padding
        self.rescue_mode = rescue_mode
        self.cache_dir = Path(cache_dir) if cache_dir is not None else None
        self.force_cache = force_cache
        self.cache_n_workers = cache_n_workers
        self._dock_fn = _dock_fn
        self._cluster_kwargs = _cluster_kwargs or {}

    def run(self, resume: bool = True, _repo: Any = None) -> ScreeningResult:
        """実行メインループ。Dask LocalCluster による並列実行。_repo はテスト用依存注入。"""
        from distributed import Client, LocalCluster

        t0 = time.monotonic()
        self.log_path.parent.mkdir(parents=True, exist_ok=True)

        if not resume and self.hdf5_path.exists():
            self.hdf5_path.unlink()

        if _repo is None:
            from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
                HDF5DockingResultRepository,
            )
            _repo = HDF5DockingResultRepository(
                self.hdf5_path,
                mode="append",
                schema_version=self.schema_version,
            )
        repo = _repo

        all_pairs = list(self._enumerate_pairs())
        total = len(all_pairs)
        unprocessed = self._filter_unprocessed(repo)
        reused = total - len(unprocessed)

        pairs_by_protein: Dict[str, List[int]] = {}
        for protein_id, compound_index in unprocessed:
            pairs_by_protein.setdefault(protein_id, []).append(compound_index)

        new_pairs = 0
        failed = 0

        with open(self.log_path, "a") as log_fp:
            # grid_box missing のペアをログに記録し除外
            valid_pairs_by_protein: Dict[str, List[int]] = {}
            for protein_id, compound_indices in pairs_by_protein.items():
                protein = self.protein_set[protein_id]
                grid_box = self.grid_box_cache.get(protein)
                if grid_box is None:
                    if self.grid_box_missing_policy == "error":
                        raise ValueError(f"GridBox missing for protein '{protein_id}'")
                    for compound_index in compound_indices:
                        log_fp.write(json.dumps({
                            "protein_id": protein_id,
                            "compound_index": compound_index,
                            "status": "failed",
                            "score": None,
                            "elapsed_sec": 0.0,
                            "error": "grid_box_missing",
                        }) + "\n")
                else:
                    valid_pairs_by_protein[protein_id] = compound_indices

            # Phase 1: receptor cache pre-generation (unidock2 only)
            if (valid_pairs_by_protein
                    and self.cache_dir
                    and self.backend == 'unidock2'):
                cache_results = self.prepare_receptor_caches(
                    valid_pairs_by_protein,
                )
                for p_hash, info in cache_results.items():
                    log_fp.write(json.dumps({
                        "phase": "cache_prep",
                        "protein_content_hash": p_hash,
                        "cache_path": info.get("cache_path"),
                        "error": info.get("error"),
                    }) + "\n")

            # Phase 2: docking
            if valid_pairs_by_protein:
                cluster_kwargs = {
                    "n_workers": self.dask_n_workers,
                    "threads_per_worker": 1,
                    "memory_limit": "4GB",
                    **self._cluster_kwargs,
                }
                cluster = LocalCluster(**cluster_kwargs)
                client = Client(cluster)
                try:
                    futures = self._submit_to_dask(client, valid_pairs_by_protein)
                    n, f = self._collect_and_save(futures, repo, log_fp)
                    new_pairs += n
                    failed += f
                finally:
                    client.close()
                    cluster.close()

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

    def prepare_receptor_caches(
        self,
        pairs_by_protein: Dict[str, List[int]],
    ) -> Dict[str, dict]:
        """Dask で UniDock2 receptor cache を並列生成する。

        既存 cache がある受容体はスキップ (force_cache=True で上書き)。
        cache_dir が None または backend が unidock2 でない場合は空 dict を返す。

        Args:
            pairs_by_protein: protein_id -> compound_indices (対象タンパク質の列挙用)。

        Returns:
            {protein_content_hash: {"cache_path": str|None, "error": str|None}}
        """
        if self.cache_dir is None or self.backend != 'unidock2':
            return {}

        from distributed import Client, LocalCluster, as_completed

        self.cache_dir.mkdir(parents=True, exist_ok=True)
        protein_hashes = self.protein_set.content_hashes()

        # 対象受容体のリストを構築 (content_hash で重複排除)
        seen_hashes: Dict[str, str] = {}  # content_hash -> protein_id
        tasks: List[dict] = []
        for protein_id in pairs_by_protein:
            protein = self.protein_set[protein_id]
            p_hash = protein_hashes[protein_id]
            if p_hash in seen_hashes:
                continue
            seen_hashes[p_hash] = protein_id

            cache_path = self.cache_dir / f"{p_hash}.json"
            if cache_path.exists() and not self.force_cache:
                continue

            grid_box = self.grid_box_cache.get(protein)
            if grid_box is None:
                continue

            padded_size = [float(s) + 2.0 * self.extra_padding for s in grid_box.size]
            tasks.append({
                "protein_path": str(protein.path),
                "protein_content_hash": p_hash,
                "grid_center": [float(c) for c in grid_box.center],
                "grid_size": padded_size,
            })

        if not tasks:
            return {}

        n_workers = self.cache_n_workers or self.dask_n_workers
        cluster = LocalCluster(
            n_workers=n_workers,
            threads_per_worker=1,
            memory_limit="4GB",
            **self._cluster_kwargs,
        )
        client = Client(cluster)
        results: Dict[str, dict] = {}
        try:
            futures = {}
            for t in tasks:
                f = client.submit(
                    _prepare_cache_for_dask,
                    t["protein_path"],
                    t["protein_content_hash"],
                    t["grid_center"],
                    t["grid_size"],
                    str(self.cache_dir),
                    self.force_cache,
                )
                futures[f] = t["protein_content_hash"]

            for future in as_completed(list(futures.keys())):
                p_hash = futures[future]
                try:
                    result = future.result()
                    results[p_hash] = result
                except Exception as e:
                    results[p_hash] = {
                        "protein_content_hash": p_hash,
                        "cache_path": None,
                        "error": f"{type(e).__name__}: {e}",
                    }
        finally:
            client.close()
            cluster.close()

        return results

    def _submit_to_dask(
        self,
        client: Any,
        pairs_by_protein: Dict[str, List[int]],
    ) -> Dict[Any, str]:
        """各タンパク質の未処理ペアをDaskに submit する。

        Returns dict[future → protein_id].
        """
        protein_hashes = self.protein_set.content_hashes()
        dock_fn = self._dock_fn if self._dock_fn is not None else dock_one_protein_via_tool
        futures: Dict[Any, str] = {}

        for protein_id, compound_indices in pairs_by_protein.items():
            protein = self.protein_set[protein_id]
            grid_box = self.grid_box_cache.get(protein)
            if grid_box is None:
                raise ValueError(f"GridBox が grid_box_cache 未登録です: protein_id={protein_id}")
            p_hash = protein_hashes[protein_id]
            compound_hashes = {
                idx: self.compound_set.get_compound_hash(idx)
                for idx in compound_indices
            }

            padded_size = [float(s) + 2.0 * self.extra_padding for s in grid_box.size]
            future = client.submit(
                dock_fn,
                str(protein.path),
                protein_id,
                p_hash,
                str(self.compound_set.path),
                compound_indices,
                compound_hashes,
                [float(c) for c in grid_box.center],
                padded_size,
                exhaustiveness=self.exhaustiveness,
                top_n_poses=self.top_n_poses,
                backend=self.backend,
                search_mode=self.search_mode,
                rescue_mode=self.rescue_mode,
                cache_dir=str(self.cache_dir) if self.cache_dir else None,
            )
            futures[future] = protein_id

        return futures

    def _collect_and_save(
        self,
        futures: Dict[Any, str],
        repo: Any,
        log_fp: Any,
    ) -> tuple[int, int]:
        """futures を as_completed で受信し、結果をHDF5に保存してJSONLに記録する。

        Returns (new_pairs, failed_pairs).
        """
        from distributed import as_completed

        new_pairs = 0
        failed = 0
        use_bundle = getattr(repo, "schema_version", "v2") == "v3"
        protein_hashes = self.protein_set.content_hashes() if use_bundle else {}

        for future in as_completed(list(futures.keys())):
            protein_id = futures[future]
            try:
                results = future.result()

                if use_bundle:
                    bundle_entries = []
                    for r in results:
                        if r.get("error"):
                            failed += 1
                            log_fp.write(json.dumps({
                                "protein_id": r["protein_id"],
                                "compound_index": r["compound_index"],
                                "status": "failed",
                                "score": None,
                                "elapsed_sec": r["elapsed_sec"],
                                "error": r.get("error"),
                            }) + "\n")
                        else:
                            bundle_entries.append({
                                "compound_hash": r["compound_content_hash"],
                                "score": r["score"],
                                "pose_blob": r.get("pose_blob"),
                                "source": self.backend,
                                "top_n": self.top_n_poses,
                            })
                            new_pairs += 1
                            log_fp.write(json.dumps({
                                "protein_id": r["protein_id"],
                                "compound_index": r["compound_index"],
                                "status": "new",
                                "score": r["score"],
                                "elapsed_sec": r["elapsed_sec"],
                            }) + "\n")
                        log_fp.flush()
                    if bundle_entries:
                        p_hash = protein_hashes[protein_id]
                        repo.write_bundle(p_hash, bundle_entries)
                else:
                    for r in results:
                        if r.get("error"):
                            failed += 1
                            log_fp.write(json.dumps({
                                "protein_id": r["protein_id"],
                                "compound_index": r["compound_index"],
                                "status": "failed",
                                "score": None,
                                "elapsed_sec": r["elapsed_sec"],
                                "error": r.get("error"),
                            }) + "\n")
                        else:
                            if r.get("pose_blob") is not None:
                                self._save_to_hdf5(
                                    r["protein_content_hash"],
                                    r["protein_id"],
                                    r["compound_index"],
                                    r["compound_content_hash"],
                                    r["score"],
                                    r["pose_blob"],
                                )
                            new_pairs += 1
                            log_fp.write(json.dumps({
                                "protein_id": r["protein_id"],
                                "compound_index": r["compound_index"],
                                "status": "new",
                                "score": r["score"],
                                "elapsed_sec": r["elapsed_sec"],
                            }) + "\n")
                        log_fp.flush()

            except Exception as e:
                failed += 1
                log_fp.write(json.dumps({
                    "protein_id": protein_id,
                    "compound_index": -1,
                    "status": "failed",
                    "score": None,
                    "elapsed_sec": 0.0,
                    "error": str(e),
                }) + "\n")
                log_fp.flush()

        return new_pairs, failed

    def _enumerate_pairs(self) -> Iterator[tuple[str, int]]:
        """(protein_id, compound_index) の全組合せを yield"""
        for protein in self.protein_set:
            for i in range(self.compound_set.get_compound_count()):
                yield (protein.id, i)

    def _filter_unprocessed(self, repo: Any) -> List[tuple[str, int]]:
        """HDF5 既存キーでフィルタし、未実行ペアのみ返す。"""
        protein_hashes = self.protein_set.content_hashes()
        use_bundle = getattr(repo, "schema_version", "v2") == "v3"

        if use_bundle:
            existing_keys = repo.get_all_keys_bundle()
            result = []
            for protein_id, compound_index in self._enumerate_pairs():
                p_hash = protein_hashes[protein_id]
                c_hash = self.compound_set.get_compound_hash(compound_index)
                if (p_hash, c_hash) not in existing_keys:
                    result.append((protein_id, compound_index))
            return result

        result = []
        for protein_id, compound_index in self._enumerate_pairs():
            p_hash = protein_hashes[protein_id]
            c_hash = self.compound_set.get_compound_hash(compound_index)
            if not repo._exists(p_hash, c_hash):
                result.append((protein_id, compound_index))
        return result

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
        import numpy as np

        self.hdf5_path.parent.mkdir(parents=True, exist_ok=True)
        with h5py.File(self.hdf5_path, "a", libver="latest") as f:
            f.swmr_mode = True
            group_path = f"/results/{protein_hash}/{compound_hash}"
            if group_path in f:
                return
            g = f.require_group(group_path)
            g.attrs["protein_id"] = protein_id
            g.attrs["compound_index"] = compound_index
            dt_str = h5py.string_dtype(encoding="utf-8")
            g.create_dataset("docking_score", data=np.float32(score))
            g.create_dataset("pose_blob", data=np.frombuffer(pose_blob, dtype=np.uint8))
            g.create_dataset(
                "computed_at",
                data=datetime.now(timezone.utc).isoformat(),
                dtype=dt_str,
            )
            g.create_dataset("source", data="vina", dtype=dt_str)
            g.create_dataset("top_n", data=np.int32(self.top_n_poses))
