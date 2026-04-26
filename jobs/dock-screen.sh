#!/bin/bash
# 任意サイズの virtual screening ジョブ (unidock2 + DMS cache 利用).
#
# 入力:
#   - 受容体 DMS:  results/receptors_dms/*.dms (jobs/prep-receptors{,-array}.sh で事前生成)
#   - リガンド SDF: jobs/inputs/ligands_sdf/*.sdf
#   - グリッド:     jobs/inputs/boxes.tsv (receptor 一覧と一致)
# 出力: results/dock-screen.<JOB_ID>/<receptor_id>/poses.sdf + unidock.log
#       logs/dock-screen.<JOB_ID>.results.tar.gz (pull で取得)
#
# 並列方針: 8 並列 (bash background + wait), no MPS
# 性能 (実測): node_q × 8 並列で ~6 秒/受容体 (100 ligands batch)
#   - 100 受容体 → 10 分
#   - 1000 受容体 → ~100 分
#
# Resume: 既に poses.sdf が生成済の受容体は skip する (途中再開対応)

#$ -cwd
#$ -l node_q=1
#$ -l h_rt=03:30:00
#$ -N dock-screen
#$ -o logs/dock-screen.$JOB_ID.out
#$ -e logs/dock-screen.$JOB_ID.err

set -eu
source $HOME/.bashrc

SIF=".sif/tsubame-env.sif"
INPUTS="jobs/inputs"
DMS_DIR="results/receptors_dms"
BOXES_TSV="$INPUTS/boxes.tsv"
RESULTS="results/dock-screen.$JOB_ID"
N_PARALLEL=8

# Resume 対象を指定したい場合は環境変数で渡す (デフォルト = 新規ディレクトリ)
RESUME_FROM="${RESUME_FROM:-}"
if [[ -n "$RESUME_FROM" && -d "$RESUME_FROM" ]]; then
    RESULTS="$RESUME_FROM"
    echo "RESUMING from existing $RESULTS"
fi

if [[ ! -f "$SIF" ]]; then echo "ERROR: sif not found at $SIF" >&2; exit 1; fi
if [[ -z "$(ls "$DMS_DIR"/*.dms 2>/dev/null || true)" ]]; then
    echo "ERROR: receptor DMS files not found in $DMS_DIR" >&2
    echo "ERROR: run jobs/prep-receptors{,-array}.sh first" >&2
    exit 1
fi
mkdir -p "$RESULTS"

echo "=== job info ==="
echo "host       : $(hostname)"
echo "JOB_ID     : ${JOB_ID:-<unset>}"
echo "nvidia-smi : $(nvidia-smi -L | head -1 || echo none)"
echo "results    : $RESULTS"
echo "started    : $(date)"

LIGANDS_BATCH="$RESULTS/_ligands_batch.txt"
ls "$INPUTS/ligands_sdf"/*.sdf > "$LIGANDS_BATCH"
N_LIGANDS=$(wc -l < "$LIGANDS_BATCH")
N_RCP=$(tail -n +2 "$BOXES_TSV" | wc -l)
echo "ligands    : $N_LIGANDS 件"
echo "receptors  : $N_RCP 件 (boxes.tsv)"
echo "parallel   : $N_PARALLEL (no MPS)"

echo "=== docking 並列実行 ==="
i=0 done=0 skipped=0 missing_dms=0
while IFS=$'\t' read -r rcp_id cx cy cz sx sy sz; do
    i=$((i+1))
    rcp="$DMS_DIR/$rcp_id.dms"
    out_dir="$RESULTS/$rcp_id"
    out_sdf="$out_dir/poses.sdf"

    # Resume: 既に処理済ならスキップ
    if [[ -s "$out_sdf" ]] && grep -q '^\$\$\$\$' "$out_sdf" 2>/dev/null; then
        skipped=$((skipped+1))
        continue
    fi
    if [[ ! -f "$rcp" ]]; then
        echo "  SKIP $rcp_id: DMS missing"
        missing_dms=$((missing_dms+1))
        continue
    fi
    mkdir -p "$out_dir"

    cat > "$out_dir/config.yaml" <<EOF
target_center: [$cx, $cy, $cz]
target_size: [$sx, $sy, $sz]
output_docking_pose_sdf_file_name: $out_sdf
EOF

    (
        apptainer exec --nv -B "$PWD:/work" --pwd /work "$SIF" \
            /opt/conda/envs/unidock2/bin/unidock2 docking \
                -r "$rcp" \
                -lb "$LIGANDS_BATCH" \
                -c "$cx" "$cy" "$cz" \
                -o "$out_sdf" \
                -cf "$out_dir/config.yaml" \
                > "$out_dir/unidock.log" 2>&1
        rc=$?
        n_poses=$(grep -c '^\$\$\$\$' "$out_sdf" 2>/dev/null || echo 0)
        echo "  [$(date +%T)] [$i/$N_RCP] $rcp_id: rc=$rc poses=$n_poses"
    ) &

    while [[ $(jobs -rp | wc -l) -ge $N_PARALLEL ]]; do
        wait -n
    done
done < <(tail -n +2 "$BOXES_TSV")
wait

echo "=== docking 完了: $(date) ==="
echo "  resumed skip: $skipped, missing dms: $missing_dms"
N_OUT=$(find "$RESULTS" -name "poses.sdf" -size +0 2>/dev/null | xargs grep -c '^\$\$\$\$' 2>/dev/null | awk -F: '{s+=$NF} END {print s}')
echo "  total poses: $N_OUT"

TARBALL="logs/dock-screen.$JOB_ID.results.tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$RESULTS")" "$(basename "$RESULTS")"
echo "  tarball: $TARBALL ($(du -h "$TARBALL" | cut -f1))"
