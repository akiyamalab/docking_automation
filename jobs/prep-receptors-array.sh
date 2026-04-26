#!/bin/bash
# 受容体 PDB → DMS の array job 並列前処理 (1000+ 件向け)
#
# 入力: jobs/inputs/receptors_pdb/AF-*.pdb
# 出力: results/receptors_dms/<rcp_id>.dms
#
# 設計:
#   - 受容体一覧をソート → タスク数 N で chunk 分割
#   - 各 array task は cpu_40 で 40 並列に protein_prep を実行
#   - 1000 件 / 20 タスク = 50 件/タスク (40 並列で ~13s + apptainer 起動 ~10s)
#
# AGE は出力ファイル名に $JOB_ID.$TASK_ID を含めること.

#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=00:30:00
#$ -N prep-receptors-array
#$ -t 1-20
#$ -o logs/prep-receptors-array.$JOB_ID.$TASK_ID.out
#$ -e logs/prep-receptors-array.$JOB_ID.$TASK_ID.err

set -eu
source $HOME/.bashrc

SIF=".sif/tsubame-env.sif"
INPUTS="jobs/inputs/receptors_pdb"
OUTDIR="results/receptors_dms"
N_TASKS=20            # = SGE_TASK_LAST (1-20)
N_PARALLEL=40

mkdir -p "$OUTDIR"

echo "=== task info ==="
echo "host: $(hostname), JOB_ID=$JOB_ID, TASK=$SGE_TASK_ID/$N_TASKS, $(date)"

# 全受容体一覧 (sorted) を取得し、自タスクの担当範囲を計算
# 競合回避のため task-local リストに書き出す (共有ファイルだと並行書込が壊れる)
ALL_LIST="$OUTDIR/_all.${JOB_ID}.${SGE_TASK_ID}.txt"
ls -1 "$INPUTS"/AF-*.pdb | sort > "$ALL_LIST"
N_TOTAL=$(wc -l < "$ALL_LIST")

# chunk size = ceil(N_TOTAL / N_TASKS)
CHUNK=$(( (N_TOTAL + N_TASKS - 1) / N_TASKS ))
START=$(( (SGE_TASK_ID - 1) * CHUNK + 1 ))
END=$(( SGE_TASK_ID * CHUNK ))
[[ $END -gt $N_TOTAL ]] && END=$N_TOTAL
echo "全受容体: $N_TOTAL, この task の担当: $START..$END (chunk size $CHUNK)"

# 担当範囲の PDB のみを抽出した task-local list
TASK_LIST="$OUTDIR/_task.${JOB_ID}.${SGE_TASK_ID}.txt"
sed -n "${START},${END}p" "$ALL_LIST" > "$TASK_LIST"

# Python orchestrator (apptainer 1 回起動)
apptainer exec -B "$PWD:/work" --pwd /work "$SIF" \
    /opt/conda/envs/unidock2/bin/python <<PYEOF
import os, subprocess, sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

OUTDIR = Path("$OUTDIR")
N = $N_PARALLEL
TASK_LIST = Path("$TASK_LIST")

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

with open(TASK_LIST) as f:
    pdbs = [Path(line.strip()) for line in f if line.strip()]
print(f"task targets: {len(pdbs)}", flush=True)

ok, fail, skip = 0, 0, 0
with ProcessPoolExecutor(max_workers=N) as ex:
    futs = {ex.submit(prep, p): p for p in pdbs}
    for i, f in enumerate(as_completed(futs), 1):
        rid, status, size = f.result()
        if status == "done":  ok += 1
        elif status == "skip": skip += 1
        else:                  fail += 1
        if i % 10 == 0 or i == len(pdbs):
            print(f"  [{i:3d}/{len(pdbs)}] {rid}: {status}", flush=True)

print(f"summary[task=$SGE_TASK_ID]: done={ok} skip={skip} fail={fail}")
PYEOF

# task list と ALL_LIST は完了後削除
rm -f "$TASK_LIST" "$ALL_LIST"

echo "=== task=$SGE_TASK_ID 完了: $(date) ==="
