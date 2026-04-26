#!/bin/bash
# 受容体 PDB → DMS の並列前処理 (TSUBAME 側でキャッシュとして保持)
#
# 入力: jobs/inputs/receptors_pdb/*.pdb
# 出力: results/receptors_dms/<rcp_id>.dms (push の --delete から exclude されるため永続)
#
# 課金: cpu_40 × 30min = 0.15 × 0.5h ≈ 5 円相当 (実時間が短ければそれ未満)

#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=00:30:00
#$ -N prep-receptors
#$ -o logs/prep-receptors.$JOB_ID.out
#$ -e logs/prep-receptors.$JOB_ID.err

set -eu
source $HOME/.bashrc

SIF=".sif/tsubame-env.sif"
INPUTS="jobs/inputs/receptors_pdb"
OUTDIR="results/receptors_dms"
N_PARALLEL=40

mkdir -p "$OUTDIR"

echo "=== job info ==="
echo "host: $(hostname), JOB_ID=$JOB_ID, $(date)"
echo "input PDBs : $(ls "$INPUTS"/*.pdb 2>/dev/null | wc -l)"
echo "outdir     : $OUTDIR"
echo "parallel   : $N_PARALLEL"

# Python orchestrator を sif 内で実行 (apptainer exec を 1 回だけ呼ぶ)
apptainer exec -B "$PWD:/work" --pwd /work "$SIF" \
    /opt/conda/envs/unidock2/bin/python <<PYEOF
import os, subprocess, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

INPUTS = Path("$INPUTS")
OUTDIR = Path("$OUTDIR")
N = $N_PARALLEL

def prep(pdb_path):
    rcp_id = pdb_path.stem
    out = OUTDIR / f"{rcp_id}.dms"
    if out.exists() and out.stat().st_size > 0:
        return (rcp_id, "skip", 0)
    cmd = [
        "/opt/conda/envs/unidock2/bin/unidock2", "protein_prep",
        "-r", str(pdb_path),
        "-o", str(out),
    ]
    log = OUTDIR / f"{rcp_id}.prep.log"
    with open(log, "w") as f:
        rc = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT).returncode
    return (rcp_id, "done" if rc == 0 else f"rc={rc}", out.stat().st_size if out.exists() else 0)

pdbs = sorted(INPUTS.glob("*.pdb"))
print(f"prep targets: {len(pdbs)}")

ok, fail, skip = 0, 0, 0
with ProcessPoolExecutor(max_workers=N) as ex:
    futs = {ex.submit(prep, p): p for p in pdbs}
    for i, f in enumerate(as_completed(futs), 1):
        rid, status, size = f.result()
        if status == "done":  ok += 1
        elif status == "skip": skip += 1
        else:                  fail += 1
        if i % 10 == 0 or i == len(pdbs):
            print(f"  [{i:3d}/{len(pdbs)}] {rid}: {status} ({size} B)", flush=True)

print(f"summary: done={ok} skip={skip} fail={fail}")
PYEOF

echo "=== prep 完了: $(date) ==="
ls -lh "$OUTDIR"/ 2>&1 | head -5
echo "  total dms files: $(ls "$OUTDIR"/*.dms 2>/dev/null | wc -l)"
echo "  total size: $(du -sh "$OUTDIR" 2>/dev/null | cut -f1)"
