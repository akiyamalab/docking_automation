#!/usr/bin/env python3
"""Uni-Dock 2 receptor JSON cache の並列生成。

Uni-Dock 2 は受容体 1 件あたり `analyze_receptor_topology` に ~5 分 (受容体サイズで
2〜15 分のレンジ) かかる。`UniDock2Docking.prepare_receptor_cache()` を受容体ごとに呼び、
24 コアマシンなら 24 並列で全受容体の prep を ~24× 高速化する。

並列化の前提 (`memory/unidock2_overhead_cache.md` 参照):
- 各 proc で `OMP_NUM_THREADS=1` 必須 (default 24 だと 2 proc で 48 threads が 24 cores を奪い合う)
- DMS I/O は O_DIRECT/mmap ではなく通常 read なので page cache が効く
- GPU kernel は使われない (analyze_receptor_topology は CPU のみ)

Usage:
    python3 scripts/prepare_unidock2_caches.py \\
        --protein-list examples/input/afdb_mouse/protein_list.json \\
        --grid-dir output/grid_cache \\
        --cache-dir output/unidock2_receptor_cache \\
        --workers 24
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Iterable, List, Optional, Tuple

SCRIPT_DIR = Path(__file__).parent.resolve()
REPO_DIR = SCRIPT_DIR.parent
sys.path.insert(0, str(REPO_DIR))

logger = logging.getLogger('prepare_unidock2_caches')


def _prepare_one(
    protein_id: str,
    protein_path: str,
    grid_center: Tuple[float, float, float],
    grid_size: Tuple[float, float, float],
    cache_dir: str,
    force: bool,
) -> Tuple[str, Optional[str], float, Optional[str]]:
    """1 受容体の cache を生成 (worker プロセスで実行される)。

    Returns:
        (protein_id, cache_path or None, elapsed_sec, error_msg or None)
    """
    # worker でも OMP=1 を徹底
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    os.environ['OPENBLAS_NUM_THREADS'] = '1'

    t0 = time.time()
    try:
        from docking_automation.docking.grid_box import GridBox
        from docking_automation.docking.unidock2_docking import UniDock2Docking
        from docking_automation.molecule.protein import Protein

        protein = Protein(Path(protein_path), id=protein_id)
        grid_box = GridBox(center=list(grid_center), size=list(grid_size))

        tool = UniDock2Docking(cache_dir=Path(cache_dir))
        out_json = tool.cache_path_for(protein)
        tool.prepare_receptor_cache(protein, grid_box, out_json, force=force)

        return (protein_id, str(out_json), time.time() - t0, None)
    except Exception as e:
        return (protein_id, None, time.time() - t0, f'{type(e).__name__}: {e}')


def load_protein_list(path: Path) -> List[dict]:
    """protein_list.json の読み込み。形式は examples/input/afdb_mouse/protein_list.json に準拠。

    期待される形式 (list of dict):
        [{"id": "Q9Z0X1", "path": "AF-Q9Z0X1-F1-model_v4.pdb",
          "grid_center": [x,y,z], "grid_size": [30,30,30]}, ...]
    """
    with path.open() as f:
        data = json.load(f)
    if isinstance(data, dict):
        data = data.get('proteins', [])
    return data


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description='Uni-Dock 2 receptor cache parallel generator')
    p.add_argument('--protein-list', type=Path, required=True,
                   help='protein_list.json (id, path, grid_center, grid_size フィールドを含む list)')
    p.add_argument('--cache-dir', type=Path, required=True,
                   help='JSON cache 保存先')
    p.add_argument('--workers', type=int, default=os.cpu_count() or 4,
                   help='並列プロセス数 (default: nproc)')
    p.add_argument('--force', action='store_true',
                   help='既存キャッシュを上書き')
    p.add_argument('--limit', type=int, default=None,
                   help='処理する receptor 数の上限 (テスト用)')
    return p.parse_args()


def main() -> int:
    logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')
    args = parse_args()
    args.cache_dir.mkdir(parents=True, exist_ok=True)

    proteins = load_protein_list(args.protein_list)
    if args.limit:
        proteins = proteins[:args.limit]
    logger.info(f'targets: {len(proteins)} receptors, workers={args.workers}')

    tasks = []
    for p in proteins:
        pid = p.get('id') or p.get('protein_id')
        path = p.get('path') or p.get('protein_path')
        center = tuple(p['grid_center'])
        size = tuple(p.get('grid_size', [30.0, 30.0, 30.0]))
        tasks.append((pid, str(path), center, size))

    t_total = time.time()
    n_ok = 0
    n_err = 0
    with ProcessPoolExecutor(max_workers=args.workers) as ex:
        futures = [
            ex.submit(_prepare_one, pid, path, c, s, str(args.cache_dir), args.force)
            for pid, path, c, s in tasks
        ]
        for fut in as_completed(futures):
            pid, cache_path, elapsed, err = fut.result()
            if err is None:
                n_ok += 1
                logger.info(f'[OK]  {pid}  {elapsed:.1f}s  -> {cache_path}')
            else:
                n_err += 1
                logger.error(f'[ERR] {pid}  {elapsed:.1f}s  {err}')

    logger.info(
        f'done: {n_ok} ok / {n_err} err, total wall={time.time() - t_total:.1f}s '
        f'(avg per receptor with {args.workers}-parallel = '
        f'{(time.time() - t_total) / max(len(tasks), 1):.1f}s)'
    )
    return 0 if n_err == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
