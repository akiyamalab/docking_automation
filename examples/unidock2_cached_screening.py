"""Uni-Dock 2 キャッシュ運用 E2E: 受容体 cache を並列生成 → cached docking → HDF5 保存。

本スクリプトは小規模 PoC 用 (受容体 3 × ligand 10)。Phase 4 規模の運用では
`scripts/prepare_unidock2_caches.py` を単独で実行してから本 screening を別プロセスで回す。

Usage:
    cd /workspaces/20260422_mouse_docking/docking_automation
    # (事前に conda env `unidock2` を activate)
    python3 examples/unidock2_cached_screening.py --dry-run
    python3 examples/unidock2_cached_screening.py
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import List, Optional, Tuple

SCRIPT_DIR = Path(__file__).parent.resolve()
REPO_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_DIR))

logger = logging.getLogger('unidock2_cached_screening')


def _worker_dock(
    cache_json: str,
    ligand_sdf_list: List[str],
    grid_center: Tuple[float, float, float],
    grid_size: Tuple[float, float, float],
    protein_content_hash: str,
    compound_content_hashes: List[str],
    timeout_sec: float = 600.0,
    max_retries: int = 2,
) -> Tuple[str, List[dict], float, Optional[str]]:
    """worker: cached docking を 1 receptor 分実行。結果を dict のリストにして親に返す。

    `dock_with_cache_robust` 経由で timeout + retry 付き (稀な内部デッドロックに対策)。

    Returns:
        (protein_content_hash, result_dicts, elapsed, error)
    """
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    t0 = time.time()
    try:
        from docking_automation.docking.grid_box import GridBox
        from docking_automation.docking.unidock2_docking import UniDock2Docking

        grid_box = GridBox(center=list(grid_center), size=list(grid_size))
        tool = UniDock2Docking()
        results = tool.dock_with_cache_robust(
            cache_json=Path(cache_json),
            ligand_sdf_list=[Path(p) for p in ligand_sdf_list],
            grid_box=grid_box,
            protein_content_hash=protein_content_hash,
            compound_content_hashes=compound_content_hashes,
            timeout_sec=timeout_sec,
            max_retries=max_retries,
        )

        # DockingResult は HDF5 保存前に SDF 本体を読み出しておく
        result_dicts = []
        for r in results:
            sdf_content = r.result_path.read_text() if r.result_path and r.result_path.exists() else ''
            result_dicts.append({
                'score': r.docking_score,
                'sdf_content': sdf_content,
                'protein_content_hash': r.protein_content_hash,
                'compound_content_hash': r.compound_content_hash,
                'compoundset_content_hash': r.compoundset_content_hash,
                'compound_index': r.compound_index,
                'compound_set_id': r.compound_set_id,
                'protein_id': r.protein_id,
                'metadata': r.metadata,
            })

        return (protein_content_hash, result_dicts, time.time() - t0, None)
    except Exception as e:
        return (protein_content_hash, [], time.time() - t0, f'{type(e).__name__}: {e}')


def save_to_hdf5(result_dicts: List[dict], hdf5_dir: Path) -> None:
    """main プロセスで HDF5 repo に追記 (並列競合回避)。"""
    from docking_automation.docking.docking_result import DockingResult
    from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
        DockingResultRepositoryFactory,
        RepositoryType,
    )

    hdf5_dir.mkdir(parents=True, exist_ok=True)
    repo = DockingResultRepositoryFactory.create(
        RepositoryType.HDF5,
        base_directory=hdf5_dir,
        config={'mode': 'append'},
    )

    import tempfile
    for rd in result_dicts:
        if not rd.get('sdf_content'):
            continue
        tmp = Path(tempfile.mkstemp(suffix='.sdf')[1])
        tmp.write_text(rd['sdf_content'])
        result = DockingResult(
            result_path=tmp,
            protein_id=rd['protein_id'],
            compound_set_id=rd['compound_set_id'],
            compound_index=rd['compound_index'],
            docking_score=rd['score'],
            protein_content_hash=rd['protein_content_hash'],
            compound_content_hash=rd['compound_content_hash'],
            compoundset_content_hash=rd['compoundset_content_hash'],
            metadata=rd['metadata'],
        )
        repo.save(result)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser()
    ap.add_argument('--cache-dir', type=Path,
                    default=SCRIPT_DIR / 'output' / 'ud2_receptor_cache')
    ap.add_argument('--hdf5-dir', type=Path,
                    default=SCRIPT_DIR / 'output' / 'ud2_screening_hdf5')
    ap.add_argument('--workers', type=int, default=os.cpu_count() or 4)
    ap.add_argument('--dry-run', action='store_true')
    ap.add_argument('--ligand-sdf-dir', type=Path, default=None,
                    help='SDF (3D化済) 群のディレクトリ。未指定時はダミー 10 件を生成')
    return ap.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    args = parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    # ------ 1. 受容体 cache の並列生成 ------
    #    Phase 4 では scripts/prepare_unidock2_caches.py を単独実行するが、
    #    ここでは PoC のため同一スクリプト内で一括実行。
    from docking_automation.docking.grid_box import GridBox
    from docking_automation.docking.unidock2_docking import UniDock2Docking
    from docking_automation.molecule.protein import Protein

    # 既存の PDB 3 件を利用 (examples/input/alphafold/)
    pdb_dir = SCRIPT_DIR / 'input' / 'alphafold'
    proteins = [
        (p, [0.0, 0.0, 0.0])  # center は重心で代用 (本運用では fpocket 等で決定)
        for p in sorted(pdb_dir.glob('AF-*.pdb'))[:3]
    ]

    logger.info(f'preparing caches for {len(proteins)} receptors, workers={args.workers}')

    from scripts.prepare_unidock2_caches import _prepare_one

    if args.dry_run:
        logger.info('[dry-run] skipping cache prep + docking')
        return 0

    t_cache_start = time.time()
    cache_paths: List[Tuple[Protein, Path]] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(
                _prepare_one,
                p.stem,
                str(p),
                tuple(center),
                (30.0, 30.0, 30.0),
                str(args.cache_dir),
                False,
            ): p
            for p, center in proteins
        }
        for fut in as_completed(futures):
            pid, cache_path, elapsed, err = fut.result()
            if err is None:
                logger.info(f'  [cache OK]  {pid}  {elapsed:.1f}s')
                cache_paths.append((Protein(Path(str(futures[fut])), id=pid), Path(cache_path)))
            else:
                logger.error(f'  [cache ERR] {pid}  {err}')
    logger.info(f'cache generation: {time.time() - t_cache_start:.1f}s total')

    # ------ 2. ligand 準備 (3D 化済 SDF) ------
    if args.ligand_sdf_dir is None:
        logger.info('no --ligand-sdf-dir specified. generating 10 dummy ligands.')
        ligand_dir = args.cache_dir.parent / 'dummy_ligands'
        ligand_sdf_list = _generate_dummy_ligands(ligand_dir, n=10)
    else:
        ligand_sdf_list = sorted(args.ligand_sdf_dir.glob('*.sdf'))

    logger.info(f'ligand batch: {len(ligand_sdf_list)} SDF')

    # ------ 3. cached docking を receptor 毎に並列実行 → HDF5 集約 ------
    t_dock_start = time.time()
    all_result_dicts: List[dict] = []
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = {
            ex.submit(
                _worker_dock,
                str(cache_path),
                [str(p) for p in ligand_sdf_list],
                (0.0, 0.0, 0.0),
                (30.0, 30.0, 30.0),
                protein.content_hash,
                [p.stem for p in ligand_sdf_list],
            ): protein.id
            for protein, cache_path in cache_paths
        }
        for fut in as_completed(futures):
            pid = futures[fut]
            hash_, rds, elapsed, err = fut.result()
            if err is None:
                logger.info(f'  [dock OK]  {pid}  {len(rds)} results  {elapsed:.1f}s')
                all_result_dicts.extend(rds)
            else:
                logger.error(f'  [dock ERR] {pid}  {err}')
    logger.info(f'docking (cached): {time.time() - t_dock_start:.1f}s total')

    # ------ 4. HDF5 に集約保存 (メインプロセス単独書き手) ------
    save_to_hdf5(all_result_dicts, args.hdf5_dir)
    logger.info(f'saved to {args.hdf5_dir}')
    return 0


def _generate_dummy_ligands(out_dir: Path, n: int) -> List[Path]:
    """テスト用 n 件の SDF (3D 化済み) を生成。"""
    from rdkit import Chem
    from rdkit.Chem import AllChem

    out_dir.mkdir(parents=True, exist_ok=True)
    smiles = [
        'c1ccccc1', 'CCOCC', 'NC(=O)c1ccncc1', 'Oc1ccc(O)cc1', 'CN(C)c1ccc(N)cc1',
        'CC(=O)Nc1ccc(O)cc1', 'OC(=O)Cc1ccccc1', 'CCN(CC)CC', 'Nc1ccccc1C(=O)O',
        'Oc1ccc2ccccc2c1',
    ][:n]
    paths = []
    for i, smi in enumerate(smiles):
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        m = Chem.AddHs(m)
        AllChem.EmbedMolecule(m, randomSeed=42 + i)
        AllChem.MMFFOptimizeMolecule(m)
        p = out_dir / f'lig_{i:03d}.sdf'
        w = Chem.SDWriter(str(p))
        m.SetProp('_Name', f'lig_{i:03d}')
        w.write(m)
        w.close()
        paths.append(p)
    return paths


if __name__ == '__main__':
    sys.exit(main())
