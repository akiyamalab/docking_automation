#!/bin/bash
#$ -cwd
#$ -l cpu_4=1
#$ -l h_rt=00:45:00
#$ -N build-apptainer
#$ -o logs/build-apptainer.$JOB_ID.out
#$ -e logs/build-apptainer.$JOB_ID.err

# TSUBAME 上で Apptainer の .sif をビルドするジョブ。
# 入力: jobs/tsubame-env.def
# 出力: .sif/tsubame-env.sif (workdir 配下の隠しディレクトリ、/gs/bs に置かれる)
# 注: $HOME (個人 25GB quota) が枯渇するため group disk (/gs/bs) 側に配置。
#     .sif/ は push の --delete から exclude 済 (bin/tsubame 参照)

set -eu
source $HOME/.bashrc

DEF_FILE="jobs/tsubame-env.def"
OUT_DIR=".sif"
OUT_SIF="$OUT_DIR/tsubame-env.sif"

mkdir -p "$OUT_DIR"

# ビルド中の OCI レイヤキャッシュ等は $T4TMPDIR (ジョブローカル SSD, 自動削除) に置く。
# $HOME 上で展開すると共有 GPFS への大量小ファイル I/O で遅い。
export APPTAINER_CACHEDIR="${T4TMPDIR}/apptainer-cache"
export APPTAINER_TMPDIR="${T4TMPDIR}/apptainer-tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

echo "host: $(hostname)"
echo "apptainer: $(apptainer --version)"
echo "def file: $DEF_FILE"
echo "out: $OUT_SIF"

# --force: 既存 .sif を上書き
# 計算ノードでは fakeroot が使えるはず (公式 freesoft ドキュメント記載)
apptainer build --force "$OUT_SIF" "$DEF_FILE"

echo "=== build done ==="
ls -lh "$OUT_SIF"

# ジョブ出力ログから .sif の所在が分かるようにするため、
# logs/ にも sha256 と path を記録しておく
sha256sum "$OUT_SIF" > "logs/build-apptainer.$JOB_ID.sha256"
echo "$OUT_SIF" >> "logs/build-apptainer.$JOB_ID.sha256"
