"""UniDock2Docking.dock_with_cache の subprocess worker (timeout/retry ラッパー用)。

独立 Python プロセスで実行することで、内部 pathos プールのデッドロック時に
process group kill で確実にクリーンアップできる。

Invocation:
    python -m docking_automation.docking._unidock2_worker <args.json>

args.json schema:
    {
      "cache_json": "...",
      "ligand_sdf_list": ["..."],
      "grid_center": [x, y, z],
      "grid_size": [sx, sy, sz],
      "working_dir": "...",
      "docking_pose_sdf": "..."
    }
"""
from __future__ import annotations
import json
import os
import sys
from pathlib import Path


def main() -> int:
    os.environ.setdefault('OMP_NUM_THREADS', '1')
    os.environ.setdefault('MKL_NUM_THREADS', '1')
    os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')

    if len(sys.argv) < 2:
        print('usage: python -m docking_automation.docking._unidock2_worker <args.json>', file=sys.stderr)
        return 2

    with open(sys.argv[1]) as f:
        args = json.load(f)

    from unidock_processing.unidocktools.unidock_protocol_runner import (
        UnidockProtocolRunner,
    )

    working_dir = Path(args['working_dir'])
    working_dir.mkdir(parents=True, exist_ok=True)

    runner = UnidockProtocolRunner(
        receptor_file_name=args['cache_json'],
        ligand_sdf_file_name_list=args['ligand_sdf_list'],
        target_center=tuple(args['grid_center']),
        working_dir_name=str(working_dir),
        docking_pose_sdf_file_name=args['docking_pose_sdf'],
    )
    runner.run_unidock_protocol()
    return 0


if __name__ == '__main__':
    sys.exit(main())
