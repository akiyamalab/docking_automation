# Docking Automation Framework

マウスプロテオーム規模 × 化合物ライブラリの N×M バーチャルスクリーニング基盤

## 概要

Docking Automation Framework は、**マウスプロテオーム規模 (AFDB ~2×10⁴ proteins) × 化合物ライブラリ (~10⁴) の
N×M バーチャルスクリーニング基盤**です。段階的拡張、冪等性、content_hash ベース再開可能性を特徴とし、
AutoDock Vina (CPU) および Uni-Dock (GPU) を統一インターフェースで利用できます。

### 主な特徴

- **マウスプロテオーム対応**: AFDB 21,452 件の PDB 構造を UniProt ID シャーディングで管理
- **content_hash ベース冪等性**: 同一タンパク質・化合物ペアは HDF5 から再利用、差分のみ計算
- **段階的スクリーニング**: Phase 0→2→3→4 と規模を拡張しながら検証済み実測値を保持
- **resume 対応**: 途中中断からの再開可能（resume=True）
- **backend switchable**: Vina (CPU) / Uni-Dock (GPU) を設定で切替可能
- **Dask 並列化**: LocalCluster による並列処理で Phase 2 で 3.15× speedup を達成

## インストール方法

### 基本インストール

```bash
# パッケージのインストール（すべての依存関係が自動的にインストールされます）
pip install -e .
```

### 開発用ツールのインストール

```bash
# 開発用ツール（テスト、リンター、ドキュメント生成など）を含めたインストール
pip install -e .[dev]
```

### RDKit のインストール

RDKit は conda 経由でのインストールが推奨されています：

```bash
conda install -c conda-forge rdkit
```

### Uni-Dock (GPU バックエンド) のインストール

```bash
# Uni-Dock v1.1.0 バイナリを GitHub Releases から取得
# sm_80 / sm_90 対応。sm_75 (RTX 2080 SUPER 相当) は公式バイナリ非対応のためソースビルドが必要
# https://github.com/dptech-corp/Uni-Dock/releases

# libcuda symlink の確認（ホスト kernel driver 版と user-space lib の整合が必要）
ls -la /usr/lib/x86_64-linux-gnu/libcuda.so*

# 注意: cuda-compat-12-4 は GeForce 非対応のため削除推奨 (cmd_010 で対応済)
# sudo apt-get remove cuda-compat-12-4
```

### devcontainer 起動手順

```bash
# GPU を有効化して起動
docker run --gpus all ...

# libcuda symlink の整合を確認してから起動すること
```

### 追加依存パッケージ

- **dimorphite-dl 2.0.2**: `compound_pipeline.preprocess` でプロトン化状態生成に使用
- **dssp / mkdssp**: `protein_segmentation` でタンパク質二次構造解析に使用
- **h5py**: HDF5 結果リポジトリ
- **dask / distributed**: LocalCluster 並列実行

```bash
pip install dimorphite-dl h5py dask distributed
apt install dssp  # or conda install -c conda-forge dssp
```

## 主要コンポーネント

| クラス / モジュール | 役割 |
|---|---|
| `Protein` | タンパク質構造を表現。`Protein(path)` コンストラクタで直接生成 |
| `ProteinSet` | 複数タンパク質の集約。AFDB 21,452 件に対応 |
| `CompoundSet` | 複数化合物の集約 |
| `GridBox` | ドッキング探索空間の定義 |
| `GridBoxCache` | JSON 永続キャッシュ。`missing_policy` 3種: skip / error / fallback_centroid |
| `compound_pipeline.preprocess` | Dimorphite-DL 2.0.2 + RDKit ETKDGv3 (seed=42) による前処理 |
| `AutoDockVina` | AutoDock Vina を使用した CPU ドッキング |
| `UniDockDocking` | Uni-Dock GPU バックエンド。backend switchable |
| `ScreeningRunner` | N×M スクリーニング司令塔。backend=vina/unidock, rescue_mode, Dask LocalCluster, extra_padding |
| `HDF5DockingResultRepository` | HDF5 結果リポジトリ。schema v2 (legacy) + v3 (protein-bundle) |

### 基本的なワークフロー

1. **タンパク質と化合物の準備**:
   - タンパク質構造ファイル（PDB, MOL2 など）の読み込み
   - 化合物ファイル（SDF, MOL2 など）の読み込み

2. **グリッドボックスの設定**:
   - ドッキング計算を行う空間の定義
   - 既知のリガンド位置や活性部位情報を元に設定

3. **ドッキング計算の実行**:
   - ドッキングツールの選択と設定
   - 計算の実行

4. **結果の解析**:
   - スコアによるランキング
   - ポーズの可視化と評価

### コード例

```python
from docking_automation.molecule import Protein, CompoundSet
from docking_automation.docking import AutoDockVina, GridBox
from docking_automation.docking.autodockvina_docking import AutoDockVinaParameters

# 1. タンパク質と化合物の準備
protein = Protein("path/to/protein.pdb")
compounds = CompoundSet("path/to/compounds.sdf")

# 2. グリッドボックスの設定
# 結晶構造のリガンド位置を中心とする場合
grid_box = GridBox(center=(15.0, 23.0, 36.0), size=(20.0, 20.0, 20.0))

# 3. ドッキングパラメータの設定
params = AutoDockVinaParameters(
    exhaustiveness=8,  # 探索の徹底度
    num_modes=9,       # 出力するポーズの数
    energy_range=3.0   # 出力するポーズのエネルギー範囲
)

# 4. ドッキング計算の実行
docking_tool = AutoDockVina()
results = docking_tool.run_docking(protein, compounds, grid_box, params)

# 5. 結果の解析
top_hits = results.get_top(10)  # 上位10件の結果を取得

for i, result in enumerate(top_hits):
    print(f"{i+1}. Score: {result.docking_score}, Compound: {result.compound_id}")
    print(f"   Pose file: {result.result_path}")
```

より詳細な例は [examples/](examples/) を参照してください。

## 段階的拡張・冪等性・content_hash

### content_hash ベースの冪等性

`HDF5DockingResultRepository` は各ドッキングペアを `protein_content_hash + compound_content_hash` で識別します。
同一ペアは HDF5 から結果を再利用し、新規差分のみ計算するため、重複計算を完全に排除します。

```
同一ペアの判定: sha256(protein_pdb_bytes) + sha256(compound_mol_bytes)
→ 一致すればスキップ、不一致なら新規ドッキング
```

### 化合物ライブラリ拡張時の挙動

化合物ライブラリを 10 件 → 20 件に拡張した場合：

- 既存の 10 件に対応する全ペアは HDF5 から **再利用**（再計算不要）
- 新規 10 件のみドッキング計算を実行

### resume 対応

`ScreeningRunner` は途中中断からの再開をサポートします：

```python
runner = ScreeningRunner(backend="vina", resume=True)
runner.run(protein_set, compound_set, grid_box_cache)
# → 完了済みペアはスキップ、未完了分のみ計算継続
```

### HDF5 スキーマ

| スキーマ | 説明 |
|---|---|
| v2 (legacy) | 化合物別フラット構造 |
| v3 (protein-bundle) | タンパク質別バンドル構造。136GB 規模の全プロテオームスクリーニングに対応 |

## examples/ ガイド

| スクリプト | 目的 |
|---|---|
| `simple_docking.py` | 基本動作確認 (1 protein × 1 compound) |
| `grid_box_from_crystal_ligand.py` | 結晶構造リガンドから GridBox 生成 |
| `nxm_poc.py` | Phase 0 PoC (10×10 Vina 逐次) |
| `nxm_poc_parallel.py` | Phase 2 (10×10 Vina Dask 並列, 3.15× speedup) |
| `phase2_e2e.py` | Phase 2 E2E + resume 冪等性検証 |
| `phase3_gpu_e2e.py` | Phase 3 UniDock GPU E2E (10×100) |
| `phase3_vina_correlation.py` | Phase 3 Vina-UniDock 相関 r=0.9580 |
| `unidock_e2e_test.py` | Uni-Dock GPU 動作確認 (5 ペア小規模) |
| `dask_executor_example.py` | Dask LocalCluster 使用例 |
| `hdf5_repository_modes_example.py` | HDF5 schema 切替例 |
| `analyze_scores.py` | matplotlib 解析 (分布/ヒートマップ/top 受容体). HDF5 と Uni-Dock 生 PDBQT 両対応 |
| `render_top_poses.py` | PyMOL headless 描画 (top-K ポーズを PNG 出力). HDF5 と PDBQT 両対応 |
| `unidock2_cached_screening.py` | Uni-Dock 2 E2E (受容体 cache 並列生成 → cached docking → HDF5 保存) |

## Phase 毎の実測値

| Phase | 対象 | 規模 | 実測値 |
|---|---|---|---|
| Phase 0 | PoC (Vina CPU 逐次) | 10×10 | 12.15 s/pair, avg score -5.933 |
| Phase 2 | Vina 並列化 | 10×10 | **speedup 3.15×** (4 workers) |
| Phase 3 | Uni-Dock GPU (RTX 2080 SUPER sm_75) | 10×100 | 249.3 s/1000 ペア, 690 valid / 310 failed |
| Phase 3 | Vina vs UniDock 相関 | 10×10 | **Pearson r = 0.9580** (production path) |
| Phase 4 | HPC 想定 (H100×8) | 2×10⁸ | 外挿 **3-7 日** |

## Phase 0 PoC: N×M Docking Execution

### Purpose

Phase 0 PoC validates the core N×M docking automation pipeline with a 10-protein × 10-compound (100 pairs) dataset, verifying HDF5 idempotency and incremental compound expansion.

### Prerequisites

- Python 3.10+
- AutoDock Vina Python binding: `pip install vina`
- Open Babel 3.1.0+: `apt install openbabel` or `conda install -c conda-forge openbabel`
- RDKit 2023.09+: `conda install -c conda-forge rdkit`
- HDF5 support: `pip install h5py`

Install all dependencies:

```bash
pip install -e .
```

### Running the PoC

```bash
cd /path/to/docking_automation

# Run1: initial 100-pair docking (exhaustiveness=4, ~20 min)
python examples/run_nxm_poc.py \
    --proteins examples/proteins/ \
    --compounds examples/actives_subset.sdf \
    --output examples/output/nxm_poc_hdf5 \
    --exhaustiveness 4

# Run2: idempotency check — all 100 pairs should be skipped (~3.5 s)
python examples/run_nxm_poc.py \
    --proteins examples/proteins/ \
    --compounds examples/actives_subset.sdf \
    --output examples/output/nxm_poc_hdf5 \
    --exhaustiveness 4

# Run3: incremental expansion to 20 compounds
python examples/run_nxm_poc.py \
    --proteins examples/proteins/ \
    --compounds examples/actives_extended20.sdf \
    --output examples/output/nxm_poc_hdf5 \
    --exhaustiveness 4
```

### Expected Output

**Run1** computes 100 new docking pairs and writes results to an HDF5 store:

```
=== N×M ドッキングPoC サマリ ===
総ペア数:        100
新規ドッキング:  100
再利用:          0
総経過時間:      1219.1s
平均ドッキング時間: 12.15s/ペア
平均スコア:      -5.933
```

**Run2** skips all 100 previously computed pairs via HDF5 cache:

```
総ペア数:        100 / 新規ドッキング: 0 / 再利用: 100 / 経過時間: ~3.5s
```

**Run3** (with 20 compounds) computes only the 100 new pairs while reusing the original 100 cached results.

Full results are saved to:
- `examples/output/nxm_poc_hdf5/` — HDF5 result store
- `examples/output/nxm_poc_metrics.jsonl` — per-pair timing and scores
- `examples/output/nxm_poc_summary.txt` — run summary

See [docs/poc_report.md](docs/poc_report.md) for a detailed PoC analysis report.

## テスト実行方法

```bash
pytest                          # 全体 (187 passed / 31 skipped)
pytest tests/docking/           # ドッキングコアのみ
pytest tests/infrastructure/    # HDF5 / リポジトリ
pytest tests/molecule/          # ProteinSet / CompoundSet
```

> **注**: 31 SKIP は Phase 3/4 向けの未実装プレースホルダー（GPU テスト等）です。

## GPU 並列化・キャッシュ運用 (RTX 4090 での実測値)

### Uni-Dock v1.1.3 + CUDA MPS (N×NLIG グリッド)

単一プロセスでは RTX 4090 の SM を飽和できないため、**CUDA MPS** (`nvidia-cuda-mps-control -d`) による複数
プロセス並列が必須。MPS 無しの並列は sequential より遅化 (0.92×)。

Wall clock [s] (exhaustiveness=8, DMS 入力):

| NLIG\N | 1 | 2 | 4 | 8 | 16 |
|---|---|---|---|---|---|
| 10  | 51.6 | 55.9  | 68.6  | 77.6  | 113.1 |
| 20  | 81.7 | 97.6  | 111.8 | 141.5 | 226.6 |
| 40  | 95.5 | 106.1 | 130.1 | 193.3 | 281.6 |
| 80  | 88.1 | 106.3 | 154.7 | 220.2 | 352.5 |
| 160 | 95.6 | 149.4 | 186.6 | 298.7 | **OOM** |

- Speedup 最大: **N=16×NLIG=10 で 7.30×** (小バッチ×高並列で init オーバーヘッドを希釈)
- Throughput 最大: **N=8×NLIG=160 で 257 pairs/min** (約 15,400 pairs/h)
- **N=16×NLIG=160 は CUDA OOM** (24 GB 4090 上限)。運用目安 N×NLIG ≲ 1600
- 運用指針: **NLIG を 150〜200 に固定してから N=4〜8 を選ぶ** (NLIG 増加が init 償却で優先)

### Uni-Dock 破綻スコアと Vina CPU 再採点フォールバック

Uni-Dock GPU v1.1.3 は全ペア中 ~2% で `[−200, +125]` kcal/mol の破綻スコアを返す (`REMARK VINA RESULT`)。
原因は GPU BFGS が box 外に pose を動かした後の refine_step 不足で `FLT_MAX` がセットされる既知バグ
([Uni-Dock #144](https://github.com/dptech-corp/Uni-Dock/issues/144), 修正 PR
[#182](https://github.com/dptech-corp/Uni-Dock/pull/182) 未マージ)。

検証結果 (`examples/analyze_scores.py` + Vina 再採点):
- 破綻 13 ペアすべてで Vina CPU 再採点が `[-8.1, -4.4]` kcal/mol 正常値を返却 → **入力は健全、GPU 実装固有**
- 特定リガンド (剛直多環 + カルボキシレート) が全受容体で再現性ある破綻
  - lig_020/117 (torsions≤2, 縮合 4 環): 4/4 受容体で破綻
- canary プリスクリーン + Vina 再採点の 2 段構えを推奨

### Uni-Dock 2 (v0.6.1) — receptor JSON キャッシュで 196× 高速化

Uni-Dock 2 は docking 1 回あたり ~320 秒の前処理オーバーヘッドがあり、プロファイルで **99.8% が
`analyze_receptor_topology()`** (msys.LoadDMS + prepare_receptor_residue_mol_list) と判明。GPU kernel
本体は 10 ligand で 1.6 秒と極めて高速。

解決策: `engine_checkpoint: true` で `ud2_engine_inputs.json` を保存し、2 回目以降の docking では:

```python
UnidockProtocolRunner(
    receptor_file_name='receptor_cache.json',       # .json 拡張子で analyze_receptor_topology を bypass
    ligand_json_file_name='ligand_cache.json',      # (optional) ligand prep も bypass
    ...
).run_unidock_protocol()
```

| 条件 | wall |
|---|---:|
| キャッシュなし | 322 s |
| **フルキャッシュ** | **1.64 s (196×)** スコア完全一致 |

並列化については `OMP_NUM_THREADS=1` 設定で対称スケーリング (2 proc × 330s でそれぞれ N=1 baseline と一致)。
default OMP=24 だと複数プロセスで 48 threads × 24 cores のオーバーサブスクリプションが発生するため注意。

Phase 4 運用前提の実装:
- **`docking_automation.docking.unidock2_docking.UniDock2Docking`**: cache prep + cached docking を API 化
- **`scripts/prepare_unidock2_caches.py`**: `ProcessPoolExecutor` で複数受容体の cache を並列生成
  (OMP=1 を worker 内で自動設定、content_hash ベースで冪等)
- **`examples/unidock2_cached_screening.py`**: 受容体並列 prep → cached docking → HDF5 集約の E2E
  PoC (Phase 4 規模では `scripts/prepare_unidock2_caches.py` を事前実行)
- Python API 必須 (CLI は `ligand_json_file_name` を未対応)
- conda env `unidock2` 配下で実行 (pip 未配布, `http://quetz.dp.tech:8088/get/baymax` channel)

## Phase 4: 今後の対応

### HPC 実行計画

国内 HPC (ABCI 3.0 / TSUBAME 4.0 / Wisteria) での全プロテオームスクリーニング実行を予定。
詳細は [`context/docking_automation_hpc_plan.md`](context/docking_automation_hpc_plan.md) を参照。

### HPCExecutor 抽象

Phase 4 着手前に以下の Executor 実装を予定：

| Executor | 対象環境 |
|---|---|
| `DaskLocalExecutor` | ローカル / 小規模 |
| `DaskSlurmExecutor` | SLURM クラスタ (ABCI 等) |
| `DaskPBSExecutor` | PBS クラスタ (TSUBAME 等) |

### HDF5 v3 protein-bundle スキーマ

全プロテオームスクリーニング (2×10⁴ proteins × 10⁴ compounds) では HDF5 v3 protein-bundle スキーマを使用。
ストレージ見込み: **136 GB**。
詳細は [`context/docking_automation_hdf5_redesign.md`](context/docking_automation_hdf5_redesign.md) を参照。

### production path の整備状況 (cmd_011/012 完了)

- Uni-Dock penalty filter (スコア異常検出)
- rescue_mode (failed pair の Vina フォールバック)
- AFDB 21,452 件 PDB の UniProt ID シャーディング整理

## トラブルシューティング

### 実行時のエラー

#### メモリエラー

大規模な化合物セットを処理する際にメモリエラーが発生する場合は、バッチサイズを小さくして実行してください：

```python
results = docking_tool.run_docking(protein, compounds, grid_box, batch_size=50)
```

#### Uni-Dock が起動しない場合

```bash
# CUDA が認識されているか確認
nvidia-smi

# libcuda symlink の確認
ls -la /usr/lib/x86_64-linux-gnu/libcuda.so*

# sm_75 (RTX 2080 SUPER) は公式バイナリ非対応 → ソースビルドが必要
# sm_80 / sm_90 は公式バイナリで動作
```
