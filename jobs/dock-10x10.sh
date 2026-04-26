#!/bin/bash
# 10×10 GPU docking on TSUBAME via Apptainer + unidock2 (動作確認用最小ジョブ)
#
# 入力: jobs/inputs/{receptors_pdb,ligands_sdf}/*, jobs/inputs/boxes.tsv の先頭 10
# 出力: results/dock-10x10.<JOB_ID>/<receptor_id>/poses.sdf
#       logs/dock-10x10.<JOB_ID>.results.tar.gz

#$ -cwd
#$ -l gpu_1=1
#$ -l h_rt=00:30:00
#$ -N dock-10x10
#$ -o logs/dock-10x10.$JOB_ID.out
#$ -e logs/dock-10x10.$JOB_ID.err

set -eu
source $HOME/.bashrc

SIF=".sif/tsubame-env.sif"
INPUTS="jobs/inputs"
DMS_DIR="results/receptors_dms"
BOXES_TSV="$INPUTS/boxes.tsv"
RESULTS="results/dock-10x10.$JOB_ID"

if [[ ! -d "$DMS_DIR" ]] || [[ -z "$(ls "$DMS_DIR"/*.dms 2>/dev/null)" ]]; then
    echo "ERROR: receptor DMS files not found in $DMS_DIR" >&2
    echo "ERROR: run jobs/prep-receptors.sh first to generate them" >&2
    exit 1
fi

if [[ ! -f "$SIF" ]]; then
    echo "ERROR: sif not found at $SIF" >&2; exit 1
fi
mkdir -p "$RESULTS"

echo "=== job info ==="
echo "host       : $(hostname)"
echo "JOB_ID     : ${JOB_ID:-<unset>}"
echo "nvidia-smi : $(nvidia-smi -L | head -1 || echo none)"
echo "started    : $(date)"

# 先頭 10 リガンドを batch ファイルに
LIGANDS_BATCH="$RESULTS/_ligands_batch.txt"
ls "$INPUTS/ligands_sdf"/*.sdf | head -10 > "$LIGANDS_BATCH"
N_LIGANDS=$(wc -l < "$LIGANDS_BATCH")
echo "ligands: $N_LIGANDS 件"

echo "=== docking シリアル ==="
i=0
# 先頭 10 受容体のみ
tail -n +2 "$BOXES_TSV" | head -10 | while IFS=$'\t' read -r rcp_id cx cy cz sx sy sz; do
    i=$((i+1))
    rcp="$DMS_DIR/$rcp_id.dms"
    [[ -f "$rcp" ]] || { echo "  SKIP $rcp_id: DMS missing"; continue; }
    out_dir="$RESULTS/$rcp_id"
    mkdir -p "$out_dir"

    cat > "$out_dir/config.yaml" <<EOF
target_center: [$cx, $cy, $cz]
target_size: [$sx, $sy, $sz]
output_docking_pose_sdf_file_name: $out_dir/poses.sdf
EOF

    apptainer exec --nv -B "$PWD:/work" --pwd /work "$SIF" \
        /opt/conda/envs/unidock2/bin/unidock2 docking \
            -r "$rcp" \
            -lb "$LIGANDS_BATCH" \
            -c "$cx" "$cy" "$cz" \
            -o "$out_dir/poses.sdf" \
            -cf "$out_dir/config.yaml" \
            > "$out_dir/unidock.log" 2>&1
    rc=$?
    n_poses=$(grep -c '^\$\$\$\$' "$out_dir/poses.sdf" 2>/dev/null || echo 0)
    echo "  [$(date +%T)] [$i/10] $rcp_id: rc=$rc poses=$n_poses"
done

echo "=== docking 完了: $(date) ==="
N_OUT=$(find "$RESULTS" -name "poses.sdf" | xargs grep -c '^\$\$\$\$' 2>/dev/null | awk -F: '{s+=$NF} END {print s}')
echo "  total poses: $N_OUT (期待 100)"

TARBALL="logs/dock-10x10.$JOB_ID.results.tar.gz"
tar -czf "$TARBALL" -C "$(dirname "$RESULTS")" "$(basename "$RESULTS")"
echo "  tarball: $TARBALL"
