# NVIDIA CUDA MPS × apptainer 互換性調査 (2026-04-25)

## 依頼

> mps の調査・問題解決はかなり重要だと考えています。ほかの計算でも多用するからです。
> すみませんが、基礎的なところからじっくりと確認し、問題解決を目指していただきたいです。
> (続いて exprorer_msmd リポジトリの参照、TSUBAME 上の同様事例の参照、最小構成での切り分け等を依頼)

## 背景

100×100 docking で MPS を使った 8 並列 (`prep-receptors-array.sh` 仕込み版相当) を狙ったが、`apptainer exec` 内 unidock が hang し、原因不明のまま MPS 抜きで運用していた。これを根本切り分けする。

## 実施 (Phase 別)

| ジョブ | 内容 | 結果 |
|---|---|---|
| 7261883/7261896/7261898 | 4 テスト diag (A: serial, B: MPS, C: parallel, D: parallel+MPS) | bash 関数 export の罠で全部 rc=127 |
| 7261901 | apptainer に `--bind $CUDA_MPS_PIPE_DIRECTORY --env CUDA_MPS_PIPE_DIRECTORY=…` 追加 | hang 続く |
| 7261910/7261922 | bind 戦略 6 種比較 (A〜F) | A (no MPS) のみ通過、B〜F TIMEOUT |
| 7262067/7262158 | host 側 minimal CUDA テスト (cudatest.cu) で MPS 単独動作確認 | host = OK、apptainer 内 = hang |
| 7262136/7262147/7262148 | cuda module 切替 (12.6.0/12.4.0/12.8.0/13.1.1) | cuda 12.x module は libcudart 提供、MPS wrapper は cuda 13 のみ |
| 7262151/7262156 | unidock2 + cuda 13 sif で再テスト | hang 解消せず |
| 7262163/7262167 | minimal bind から段階的に追加 (--no-mount tmp / home, --writable-tmpfs, -B /dev, --userns, APPTAINERENV_) | 全 hang |
| 7262174/7262177 | exprorer_msmd 流 (`CUDA_MPS_PIPE_DIRECTORY=$TMPDIR/nvidia-mps`、bind 最小) | hang |
| 7262181/7262184 | container 内 ldd / nvidia-smi / mount 確認 + server.log 取得 | container は GPU 見える、bind/env 問題なし、MPS server は Status=ACTIVE まで進むが RPC コマンド受信せず |
| 7262197 | /dev/shm + 7 種フラグ最終比較 | 全 hang |

## 重要発見

1. **TSUBAME MPS wrapper は CUDA 13.0 ハードコード** (`/apps/t4/rhel9/cuda/wrapper/.bin/.nvidia-cuda-mps-control.cuda13.0`、どの cuda module を load しても同じ)
2. **conda-forge の unidock v1 は CUDA 12.6-12.x ビルドのみ** (recipe に `startswith("12")` 条件)
3. **dptech baymax channel の unidock2 0.6.1 は CUDA 13 ビルド存在** → v2 sif でバージョン揃えても hang する
4. **MPS server は client 接続を受け付けるが、`Static partitioning mode disabled` で停止し RPC コマンドが届かない**: shm 経由通信か namespace 由来の何かが阻害している

## 主要ファイル

- 設計用 def: `tsubame_skills/jobs/tsubame-env.def` (v1 → v2 へ更新、cuda 13 + unidock2)
- 検証用ジョブ: `diag-mps-*.sh` (10 種以上、最終的にすべて削除)
- 参考実装: 別 HPC の exprorer_msmd `t300k.sh`、`run_msmd.sh` (PBS Pro 系、TSUBAME 4 ではない)

## 結論

**TSUBAME 4 + apptainer 1.3.6 + MPS の組み合わせは、現時点で動作しない**。bind/env/sif version すべてを総当たりで試したが post-connection RPC 段階で hang する根本原因に到達できず。原因候補は:

1. memfd FD 渡しの mmap 失敗 (apptainer の namespace 隔離)
2. `/proc` の `hidepid=noaccess` による client→server プロセス情報参照失敗
3. CUDA 13 MPS protocol と apptainer 1.3.6 の version-specific 非互換

実用解として **MPS を諦め、parallel-no-MPS で運用** する方針を採用 (100×100 で serial 15 分 vs 8 並列 21 分、大規模だと並列恩恵あり)。

## 学びと変更

- **SKILL.md** に「MPS は apptainer から動かない」を恒久的に記録 (将来のセッションが同じ調査に時間を費やさないように)
- **CUDA_MPS_PIPE_DIRECTORY 変更** は公式は禁止と書いているが exprorer_msmd 流は変更していて daemon 起動は通る (ただし apptainer 接続は不可)
- **検証環境の整備**: 最小 cudatest.cu を nvcc で compile し host vs apptainer で挙動比較する手順は今後も使える
