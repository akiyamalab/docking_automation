"""ScreeningTool ABC — バーチャルスクリーニング共通基底。

N×M スクリーニング (多受容体 × 多 ligand) の共通パターンを表現する ABC。
受容体+grid の重い前処理 (map / topology JSON) をファイルキャッシュし、
多数 ligand を高速に回す形を基本形とする。

具体的な重い処理:
- AutoDock Vina (CPU): `compute_vina_maps` で atom-type 別 grid potential map を計算
  (receptor サイズに依存、典型的に 5〜30 秒)。`.C_H.map` などの複数ファイルで保存。
- Uni-Dock 2 (GPU): `analyze_receptor_topology` (msys.LoadDMS + residue mol リスト構築)
  で ~5 分。`ud2_engine_inputs.json` の `receptor` キーに保存。
- Uni-Dock 1 (GPU): CLI の `--write_maps` / `--maps` でバイナリ map cache 保存。

Phase 4 スケールの N×M スクリーニングでは N 受容体のキャッシュを 1 回だけ作成、
M ligand を `dock_with_cache` で高速処理する。24 core ノード上で OMP=1 並列化
すれば N 件の prep が ~N/24 × 受容体あたり処理時間で完了する。

抽象メソッド (ツールごとに実装):
- `prepare_receptor_cache(protein, grid_box, out_cache, force=False) -> Path`
- `dock_with_cache(cache, ligand_paths, grid_box, protein_content_hash, ...) -> List[DockingResult]`
- `_preprocess_compound_set(compound_set) -> PreprocessedCompoundSet`

具体メソッド (ABC で提供、全ツール共通):
- `run_docking(...)`: 受容体+ligand+grid から 1 shot docking (cache を tempdir 内に作成)
- `run_docking_with_reuse(...)`: HDF5 repo から再利用しつつ未計算ペアだけ docking
"""
from __future__ import annotations

import tempfile
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, List, Optional, Set

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.docking.grid_box import GridBox
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein

if TYPE_CHECKING:
    from docking_automation.docking.docking_parameters import SpecificDockingParametersABC
    from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
        HDF5DockingResultRepository,
    )


class ScreeningTool(ABC):
    """バーチャルスクリーニング用ドッキングツールの共通基底。

    実装クラス:
    - `AutoDockVina` (CPU): compute_vina_maps + write_maps で map cache
    - `UniDockDocking` (v1, GPU): CLI --write_maps / --maps で map cache
    - `UniDock2Docking` (v2, GPU): engine_checkpoint で receptor JSON cache
    """

    # --- 抽象: 受容体キャッシュ生成 ---
    @abstractmethod
    def prepare_receptor_cache(
        self,
        protein: Protein,
        grid_box: GridBox,
        out_cache: Path,
        force: bool = False,
    ) -> Path:
        """受容体+grid から重い前処理結果をキャッシュに保存。"""
        ...

    # --- 抽象: ligand 前処理 (ツール固有) ---
    @abstractmethod
    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        """CompoundSet をツール固有のファイル形式に前処理。"""
        ...

    # --- 抽象: キャッシュを使った docking ---
    @abstractmethod
    def dock_with_cache(
        self,
        cache: Path,
        ligand_paths: List[Path],
        grid_box: GridBox,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
    ) -> List[DockingResult]:
        """キャッシュを利用して多数 ligand をドッキング。

        ligand_paths の形式はツール依存 (Vina/v1: PDBQT, v2: SDF)。
        """
        ...

    # --- 具体: 高レベル API (ABC で共通実装) ---
    def run_docking(
        self,
        protein: Protein,
        compound_set: CompoundSet,
        grid_box: GridBox,
        additional_params: Optional['SpecificDockingParametersABC'] = None,
        compound_indices: Optional[Set[int]] = None,
    ) -> DockingResultCollection:
        """1 受容体 × 1 CompoundSet のドッキング (cache を内部で tempdir に作成)。

        N×M スクリーニングで同一受容体を多数 ligand に対し回す場合は、
        `prepare_receptor_cache` を 1 度だけ外で呼んで `dock_with_cache` を
        直接使う方がキャッシュ再利用できて効率的。
        """
        if compound_indices is not None:
            compound_set = compound_set.with_indices(compound_indices)

        prep_cs = self._preprocess_compound_set(compound_set)
        ligand_paths = [Path(p) for p in prep_cs.file_paths]
        hashes = [prep_cs.get_compound_hash(i) for i in range(len(ligand_paths))]

        with tempfile.TemporaryDirectory(prefix='screening_cache_') as tmp:
            cache_prefix = Path(tmp) / protein.content_hash
            self.prepare_receptor_cache(protein, grid_box, cache_prefix)
            results = self.dock_with_cache(
                cache=cache_prefix,
                ligand_paths=ligand_paths,
                grid_box=grid_box,
                protein_content_hash=protein.content_hash,
                compound_content_hashes=hashes,
            )

        collection = DockingResultCollection()
        collection.extend(results)
        return collection

    def run_docking_with_reuse(
        self,
        protein: Protein,
        compound_set: CompoundSet,
        grid_box: GridBox,
        additional_params: Optional['SpecificDockingParametersABC'],
        repository: 'HDF5DockingResultRepository',
        compound_indices: Optional[Set[int]] = None,
    ) -> DockingResultCollection:
        """HDF5 repo から既存結果を再利用しつつ未計算ペアだけ docking。

        content_hash ベースの冪等性保証。Phase 3/4 の resume 運用対応。
        """
        if compound_indices is not None:
            compound_set = compound_set.with_indices(compound_indices)

        prep_cs = self._preprocess_compound_set(compound_set)
        ligand_paths = [Path(p) for p in prep_cs.file_paths]
        hashes = [prep_cs.get_compound_hash(i) for i in range(len(ligand_paths))]

        # 既存分を repo から読み出し、未計算分を識別
        reused: List[DockingResult] = []
        new_idx: List[int] = []
        for i, ch in enumerate(hashes):
            if repository._exists(protein.content_hash, ch):
                try:
                    existing = repository.load_by_hashes(protein.content_hash, ch)
                    if existing is not None:
                        reused.append(existing)
                        continue
                except Exception:
                    pass
            new_idx.append(i)

        new_results: List[DockingResult] = []
        if new_idx:
            with tempfile.TemporaryDirectory(prefix='screening_cache_') as tmp:
                cache_prefix = Path(tmp) / protein.content_hash
                self.prepare_receptor_cache(protein, grid_box, cache_prefix)
                new_results = self.dock_with_cache(
                    cache=cache_prefix,
                    ligand_paths=[ligand_paths[i] for i in new_idx],
                    grid_box=grid_box,
                    protein_content_hash=protein.content_hash,
                    compound_content_hashes=[hashes[i] for i in new_idx],
                )
            for r in new_results:
                repository.save(r)

        collection = DockingResultCollection()
        collection.extend(reused)
        collection.extend(new_results)
        return collection
