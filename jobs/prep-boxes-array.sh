#!/bin/bash
# fpocket で 1000+ 受容体に対する grid box を array job 並列計算.
#
# 使用前提: sif rebuild 済 (tsubame-env.def に fpocket 追加).
# 入力: jobs/inputs/receptors_pdb/AF-*.pdb
# 出力: results/boxes_cache/<rcp_id>.json (個別 box) + 集約 boxes.tsv
#
# 注: ローカルで fpocket 実行可能な場合はローカル generate_boxes.py の方が手軽.
# このジョブは sif に fpocket 入った状態 (tsubame-env.def 更新後) でのみ動作.

#$ -cwd
#$ -l cpu_40=1
#$ -l h_rt=00:30:00
#$ -N prep-boxes-array
#$ -t 1-20
#$ -o logs/prep-boxes-array.$JOB_ID.$TASK_ID.out
#$ -e logs/prep-boxes-array.$JOB_ID.$TASK_ID.err

set -eu
source $HOME/.bashrc

SIF=".sif/tsubame-env.sif"
INPUTS="jobs/inputs/receptors_pdb"
OUTDIR="results/boxes_cache"
N_TASKS=20
N_PARALLEL=40

mkdir -p "$OUTDIR"

echo "task: $SGE_TASK_ID/$N_TASKS, host=$(hostname), $(date)"

ALL_LIST="$OUTDIR/_all.txt"
ls -1 "$INPUTS"/AF-*.pdb | sort > "$ALL_LIST"
N_TOTAL=$(wc -l < "$ALL_LIST")
CHUNK=$(( (N_TOTAL + N_TASKS - 1) / N_TASKS ))
START=$(( (SGE_TASK_ID - 1) * CHUNK + 1 ))
END=$(( SGE_TASK_ID * CHUNK ))
[[ $END -gt $N_TOTAL ]] && END=$N_TOTAL
TASK_LIST="$OUTDIR/_task.${JOB_ID}.${SGE_TASK_ID}.txt"
sed -n "${START},${END}p" "$ALL_LIST" > "$TASK_LIST"

apptainer exec -B "$PWD:/work" --pwd /work "$SIF" \
    /opt/conda/envs/unidock2/bin/python <<PYEOF
import json, subprocess, tempfile
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed

OUTDIR = Path("$OUTDIR")
TASK_LIST = Path("$TASK_LIST")
N = $N_PARALLEL

def run_fpocket(pdb_path):
    rid = pdb_path.stem
    out_json = OUTDIR / f"{rid}.json"
    if out_json.exists() and out_json.stat().st_size > 0:
        return (rid, "skip")
    with tempfile.TemporaryDirectory() as td:
        td = Path(td)
        # fpocket は同名 _out ディレクトリを作る -> tmp に link
        local_pdb = td / pdb_path.name
        local_pdb.symlink_to(pdb_path.resolve())
        rc = subprocess.run(
            ["/opt/conda/envs/unidock2/bin/fpocket", "-f", str(local_pdb)],
            capture_output=True
        ).returncode
        if rc != 0:
            return (rid, f"fp_rc={rc}")
        info = td / f"{pdb_path.stem}_out" / f"{pdb_path.stem}_info.txt"
        if not info.exists():
            return (rid, "no_info")
        # info.txt から rank 1 ポケットの中心を取得
        # フォーマット例: "Pocket 1 :\n  Score : ...\n  ..."  座標は別途 _atm.pdb から計算する必要あり
        # 簡易: pocket1_atm.pdb の重心を center, 30A 固定
        pocket_pdb = td / f"{pdb_path.stem}_out" / "pockets" / "pocket1_atm.pdb"
        if not pocket_pdb.exists():
            return (rid, "no_pocket1")
        xs, ys, zs = [], [], []
        for line in pocket_pdb.read_text().splitlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                xs.append(float(line[30:38])); ys.append(float(line[38:46])); zs.append(float(line[46:54]))
        if not xs:
            return (rid, "empty_pocket")
        cx, cy, cz = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
        # box size: pocket atoms の bounding box + padding 4A、最小 10、最大 30
        sx = min(30.0, max(10.0, max(xs)-min(xs)+8))
        sy = min(30.0, max(10.0, max(ys)-min(ys)+8))
        sz = min(30.0, max(10.0, max(zs)-min(zs)+8))
        out_json.write_text(json.dumps({"center": [cx, cy, cz], "size": [sx, sy, sz], "source": "fpocket"}))
        return (rid, "done")

with open(TASK_LIST) as f:
    pdbs = [Path(line.strip()) for line in f if line.strip()]
print(f"task targets: {len(pdbs)}", flush=True)

ok = skip = fail = 0
with ProcessPoolExecutor(max_workers=N) as ex:
    futs = {ex.submit(run_fpocket, p): p for p in pdbs}
    for i, f in enumerate(as_completed(futs), 1):
        rid, status = f.result()
        if status == "done": ok += 1
        elif status == "skip": skip += 1
        else: fail += 1
        if i % 10 == 0 or i == len(pdbs):
            print(f"  [{i}/{len(pdbs)}] {rid}: {status}", flush=True)
print(f"summary[task=$SGE_TASK_ID]: done=$ok skip=$skip fail=$fail")
PYEOF

rm -f "$TASK_LIST"
echo "task=$SGE_TASK_ID 完了: $(date)"
