"""Slurm 経由で UD2 cached docking を投入する E2E スモーク。

- dask_jobqueue SLURMCluster で worker を submit (sbatch 経由)
- 各 worker が UniDock2Docking.dock_with_cache_robust を実行
- 結果を集約

前提: Slurm (slurmctld + slurmd) + munge が起動済み, unidock2 conda env が
`/opt/miniforge/envs/unidock2` に存在、rec*_cache.json が事前生成済み。

Usage:
    /opt/miniforge/envs/unidock2/bin/python examples/unidock2_slurm_smoke.py
"""
from __future__ import annotations
import os, sys, time
from pathlib import Path

os.environ['OMP_NUM_THREADS'] = '1'

from dask_jobqueue.slurm import SLURMCluster
from dask.distributed import Client


def dock_on_worker(rec_idx: int, ligand_sdf_paths, grid_center):
    """worker 内で動く関数。UD2 conda env の python が実行する想定。"""
    import os
    os.environ['OMP_NUM_THREADS'] = '1'
    import sys
    sys.path.insert(0, '/workspaces/20260422_mouse_docking/docking_automation')
    from pathlib import Path
    from docking_automation.docking.unidock2_docking import UniDock2Docking
    from docking_automation.docking.grid_box import GridBox

    tool = UniDock2Docking()
    results = tool.dock_with_cache_robust(
        cache_json=Path(f'/tmp/ud2_test/rec{rec_idx}_cache.json'),
        ligand_sdf_list=[Path(p) for p in ligand_sdf_paths],
        grid_box=GridBox(center=list(grid_center), size=[30.0, 30.0, 30.0]),
        protein_content_hash=f'slurm_rec{rec_idx}',
        timeout_sec=180,
        max_retries=1,
    )
    return [(r.compound_index, r.docking_score) for r in results]


CENTERS = {
    0: [-2.0023, -1.3722, -1.7703],
    1: [-3.9966, 6.5890, -6.0440],
    2: [-21.8810, 1.8996, 2.2982],
}

# UD2 env の python を使うようにして Slurm worker を起動
cluster = SLURMCluster(
    queue='local',
    cores=2,
    memory='4GB',
    walltime='00:10:00',
    python='/opt/miniforge/envs/unidock2/bin/python',
    job_extra_directives=['--nodes=1'],
    job_script_prologue=['export CUDA_MPS_PIPE_DIRECTORY=/tmp/nvidia-mps',
                         'export OMP_NUM_THREADS=1'],
)
cluster.scale(jobs=3)
print('cluster:', cluster)

client = Client(cluster)
print('waiting for 3 workers via SLURM...')
client.wait_for_workers(3, timeout=120)
print(f'workers ready: {len(client.ncores())}')

ligand_paths = sorted(Path('/tmp/ud2_test/ligands_sdf').glob('lig_*.sdf'))[:10]
lig_str = [str(p) for p in ligand_paths]

t0 = time.time()
futures = [client.submit(dock_on_worker, rec_idx, lig_str, CENTERS[rec_idx])
           for rec_idx in (0, 1, 2)]
results = client.gather(futures)
wall = time.time() - t0

print(f'\n=== RESULTS (wall {wall:.1f}s via Slurm) ===')
for rec_idx, res in enumerate(results):
    print(f'rec{rec_idx}: {len(res)} ligands docked')
    for lig_idx, score in res[:3]:
        print(f'  lig_{lig_idx:03d}: score={score:.2f}')

client.close()
cluster.close()
print('DONE')
