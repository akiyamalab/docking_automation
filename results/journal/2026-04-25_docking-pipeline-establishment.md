# TSUBAME ドッキングパイプライン構築 (2026-04-25)

## 依頼

> tsubame_skills を利用して、安全に、TSUBAME 上で 10×10 などの小規模実験をおこなってください。
> (途中経過で 100×100 への拡大、unidock2 + cuda 13 への移行も指示)

## 実施

| ジョブ ID | スクリプト | リソース | 内容 | 結果 |
|---|---|---|---|---|
| 7260506 | example_cpu_small | cpu_4 trial | 疎通確認 | OK |
| 7260583 | build-apptainer | cpu_4 × tga-pharma | 初版 sif (rdkit のみ) | 1.2GB |
| 7260669 | build-apptainer (vina+) | cpu_4 | vina pip 依存で boost 不足 → 失敗 | def 修正 |
| 7260673 | build-apptainer | cpu_4 | unidock v1 (cuda 12) sif | OK 1.2GB |
| 7260678 | dock-10x10 | gpu_1 trial | 10×10 v1 docking | 100 poses ✅ |
| 7260682 | dock-10x10 (tarball) | gpu_1 trial | tarball 化追加 | OK |
| 7261126 | dock-100x100 (node_o + MPS) | node_o | MPS daemon 起動失敗 (permission denied) | ✗ |
| 7261389 | dock-100x100 (node_q + MPS) | node_q | MPS hang | ✗ |
| 7261436 | dock-100x100 (no MPS) | node_q | 並列 no-MPS、途中 cancel | — |
| 7261518 | dock-100x100 (proper MPS) | node_q | `module load cuda` 後でも hang | ✗ |
| 7261784 | dock-100x100 (serial) | gpu_1 | v1 serial | **15 min, 10000 poses** |
| 7261903 | dock-100x100 (parallel no-MPS) | node_q | v1 8 並列 | 21 min (GPU 競合) |
| 7262156 | build-apptainer (v2) | cpu_4 | unidock2 0.6.1 + cuda 13 + rdkit + matplotlib | 1.4GB ✅ |
| 7262286 | dock-10x10 (v2) | gpu_1 | v2 動作確認 | 1000 poses ✅ |
| 7262438 | prep-receptors | cpu_40 | 100 受容体 PDB→DMS 変換 | 13 秒 ✅ |
| 7262439 | dock-100x100 (v2 + DMS) | node_q | v2 + DMS cache 利用 | **10 min, 99988 poses** |

## 主要ファイル

- 入力: `docking_automation/examples/input/afdb_mouse/protein_list.json` (100 件), `actives_subset100.sdf`
- 設定: `tsubame_skills/jobs/{tsubame-env.def, dock-10x10.sh, dock-100x100.sh, prep-receptors.sh, prepare_inputs.py}`
- sif: `~/apptainer/tsubame-env.sif` (この時点では home 上、後日 `.sif/` へ移動)
- 結果: `dock-100x100.7261784.results.tar.gz` (29 MB), `dock-100x100.7262439.results.tar.gz` (33 MB)

## 学びと変更

- **trial mode は h_rt ≤ 3 分制約**: build-apptainer や docking ジョブには使えず、`tga-pharma` group 必須
- **conda-forge unidock v1 は cuda 12.x build のみ**: TSUBAME MPS daemon (cuda 13) と protocol 不整合
- **解決策**: dptech の baymax channel に **unidock2 の cuda 13 build がある** → v2 へ migration
- **DMS キャッシュは決定的**: v2 で `unidock2 protein_prep` を 1 度実行しておけば dock-screen で 30s/受容体 → 6s/受容体 に短縮
- **並列化の落とし穴**: GPU 1 枚に 8 process MPS なしで載せると競合してシリアルより遅い (15→21 min)

## 結論

unidock2 + cuda 13 + DMS キャッシュ + シリアル/並列 (MPS なし) 構成で、100×100 を **10 分で完走** できる実用基盤が整った。MPS 利用は別途検証 → [2026-04-25_mps-investigation.md](2026-04-25_mps-investigation.md) 参照。
