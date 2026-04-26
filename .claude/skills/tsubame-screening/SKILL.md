---
name: tsubame-screening
description: docking_automation を TSUBAME4 上で大規模に走らせるための docking 固有ワークフロー (受容体 PDB→DMS prep、fpocket 必須、unidock2 + cuda 13 sif、dock-screen.sh の resume 機構)。発火: docking, virtual screening, unidock2, DMS cache, fpocket, dock-screen, 大規模 docking, prep-receptors, 受容体 prep, 1000x100, screening。前提として tsubame skill (汎用) と tsubame-apptainer skill を併用。
---

# tsubame-screening skill

`docking_automation` の virtual screening を TSUBAME4 上で大規模に流すときの **docking 固有の運用知見**。
TSUBAME 操作の汎用部分 (`tsubame submit`, AGE 構文, MPS など) は `tsubame` skill 参照。
sif build は `tsubame-apptainer` skill 参照。

## いつ使うか

- 「TSUBAME で docking 流して」「大規模 screening」「N×M docking」
- 「unidock2 で hang した」「DMS cache どこ?」「fpocket box 作って」
- `docking_automation/jobs/dock-screen.sh` 等を編集・運用するとき

## ワークフロー全体像

```
[ローカル]
  1. select_1000_receptors.py   data/afdb/pdb/ から N 件選定 → inputs/receptors_pdb/
  2. prepare_inputs.py          ligand SDF を 1 分子ずつ分割 → inputs/ligands_sdf/
  3. generate_boxes.py (★fpocket)  inputs/boxes.tsv を生成
  4. tsubame push               一式を TSUBAME workdir に同期

[TSUBAME]
  5. prep-receptors-array.sh    PDB → DMS (受容体キャッシュ、永続)
  6. dock-screen.sh             N×M docking 本実行 (resume 可能)
  7. (失敗時) RESUME_FROM=results/dock-screen.<JOB_ID> dock-screen.sh で続きから

[ローカル]
  8. tsubame pull               logs/ + tarball を docking_automation/results/jobs/<JOBID>/ へ取得
```

## 重要な docking 固有の知見

### ★ Grid box は fpocket で計算する (centroid box は罠)

`unidock2` は受容体表面の空白領域に grid box が置かれると **無限探索ループに入って hang** する。
1000×100 ドッキングで 55 受容体ハングした事案あり (2026-04 検証、journal 参照)。

- **正解**: fpocket で rank 1 ポケットを推定して box にする → `generate_boxes.py` がこれを実装
- **罠**: CA 重心 + 30Å 固定 box (centroid fallback)。簡単だが大規模では一定割合 hang する

ローカルに fpocket がない場合は conda-forge から `mamba install -c conda-forge fpocket` で入る。
あるいは sif (tsubame-env.sif) に同梱して TSUBAME 上で並列計算する `prep-boxes-array.sh` 雛形もあり (sif rebuild 必要)。

### ★ 受容体 DMS prep はキャッシュする

`unidock2 docking -r receptor.pdb` を直接呼ぶと毎回 protein_prep が走り **30-50 秒/受容体**。
`unidock2 protein_prep -r X.pdb -o X.dms` で **DMS 形式に変換**しておけば docking 時は **5-10 秒/受容体**になる。

- 配置: TSUBAME workdir 上 `results/receptors_dms/<receptor_id>.dms` (push の `--delete` 対象外、永続)
- 並列実行: `prep-receptors-array.sh` (cpu_40 × 20 array、約 1 分で 1000 件完走)

### ★ unidock v1 vs unidock2 (v2)

| | unidock v1 | unidock2 |
|---|---|---|
| 入手元 | conda-forge | dptech baymax (`http://quetz.dp.tech:8088/get/baymax`) |
| CUDA | 12.x のみ | **13.x build あり** (TSUBAME MPS と整合) |
| CLI | `unidock --receptor X --gpu_batch ...` | `unidock2 docking -r X -lb LIG_LIST -c ...` |
| 入力 | PDBQT (要 meeko 変換) | PDB or DMS、SDF をそのまま |
| 推奨 | レガシー | **本リポは v2 採用** |

sif は v2 + cuda 13 で build (`docking_automation/jobs/tsubame-env.def`)。

### ★ dock-screen.sh の resume 機構

`RESUME_FROM=results/dock-screen.<前回の_JOB_ID>` を環境変数で渡すと、同じディレクトリに不足受容体だけ追記する。
完走済の `poses.sdf` (末尾 `$$$$` あり) は skip する。

`tsubame submit` ラッパーは env 渡しに対応していないので、現状は dock-screen.sh の `RESUME_FROM=` 行を一時的に
ハードコードして再投入する運用 (作業後に空に戻すこと)。

## 主要ファイル

### ローカル (`docking_automation/jobs/`)
- `prepare_inputs.py` — リガンド SDF 分割 + 受容体 PDB コピー
- `select_1000_receptors.py` — afdb tree から N 件選定 (size フィルタ付)
- `generate_boxes.py` — **fpocket で boxes.tsv 生成** (これが docking 安定化の鍵)
- `prep-receptors-array.sh` — DMS prep array job
- `dock-screen.sh` — N×M docking 本体 (resume 機構付き、node_q 8 並列、no MPS)
- `dock-10x10.sh` — smoke test 用 (gpu_1 シリアル)
- `tsubame-env.def` — sif build 定義 (unidock2 + cuda 13 + rdkit + matplotlib)

### TSUBAME workdir 上の永続データ
- `.sif/tsubame-env.sif` — 1.4 GB sif
- `results/receptors_dms/*.dms` — DMS cache (push の delete 対象外)

## 既存ジャーナル参照

実例 (`docking_automation/results/journal/`):
- `2026-04-25_docking-pipeline-establishment.md` — 10×10 から 100×100 までの基盤確立
- `2026-04-25_mps-investigation.md` — MPS+apptainer は不可、no-MPS で運用
- `2026-04-26_1000x100-virtual-screen.md` — 1000×100 完走 (fpocket box, 999,509 poses)

## エラー対応

| 症状 | 原因 | 対処 |
|---|---|---|
| `dock-screen` の一部受容体が hang | centroid box が空白領域 | `generate_boxes.py` で fpocket box 再生成 + resume |
| `dock-screen` 「DMS missing」スキップ多数 | `prep-receptors-array.sh` が一部取りこぼし | 同 array job を再投入 (resume 機能あり、`$ALL_LIST` の競合は task-local list で解決済) |
| `tsubame pull` 後の tarball が小さい | AGE が job を qstat から消すタイミングと tar 完了に race | 数秒待って `tsubame pull` 再実行 |
| `unidock2 docking` が即終了で poses 空 | 入力 SDF が無効 / DMS と box が不整合 | `unidock.log` で具体エラー確認 |

## 関連スキル
- `tsubame` — TSUBAME 共通操作 (push/pull/submit、AGE 構文、MPS、resume パターン汎用版)
- `tsubame-apptainer` — sif build
