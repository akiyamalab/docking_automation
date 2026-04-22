# Uni-Dock 統合調査レポート

**作成日**: 2026-04-22
**対象**: docking_automation Phase 3 (cmd_010 subtask_010_gunshi Part C)
**作成**: 軍師

---

## 1. 背景

Phase 4 本番規模（マウスプロテオーム 2×10⁴ × 化合物 10⁴ = 2×10⁸ ペア）を目標とする VS 基盤において、AutoDock Vina CPU 単体では計算時間が現実的でない。**Uni-Dock（DP Technology）**は AutoDock Vina 互換スコア関数を GPU で並列実行し、巨大バッチ処理を可能にする OSS ドッキングツール。本書は Phase 3 着手前の事前調査として、CLI 仕様・既存 `DockingToolABC` との整合性・統合方針を整理する。

本調査は公開情報（GitHub / arXiv 論文 / pip パッケージメタデータ）に基づく設計検討であり、実機動作検証は Phase 3 実装タスクで実施する。

---

## 2. Uni-Dock CLI 仕様

### リポジトリ
- GitHub: `dptech-corp/Uni-Dock`（公開 OSS、Apache-2.0）
- 論文: *Yu et al., "Uni-Dock: GPU-Accelerated Docking Enables Ultralarge Virtual Screening", JCTC 19, 3336-3345 (2023)*
- v1.1 系が安定、v2.0 系は開発中。Phase 3 は v1.1 を前提とする。

### 入出力形式
- **入力（受容体）**: PDBQT（AutoDock Vina と同じ）
- **入力（リガンド）**: PDBQT（複数バッチ可、`--ligand_index` ファイルで列挙）
- **出力（ポーズ）**: PDBQT（リガンドごとに 1 ファイル、`--dir` 配下）
- **スコア出力**: 標準出力 or `--score_file`（CSV 可）
- **既存 AutoDock Vina との互換性**: 同一 PDBQT スキーマで相互運用可能。MoleculeConverter の `protein_to_pdbqt` / `compound_to_pdbqt` をそのまま利用できる。

### 典型的な CLI 呼び出し（バッチモード）

```bash
unidock \
  --receptor         protein.pdbqt \
  --ligand_index     ligands.txt \
  --center_x 10.0 --center_y 10.0 --center_z 10.0 \
  --size_x   20.0 --size_y   20.0 --size_z   20.0 \
  --exhaustiveness   1 \
  --num_modes        3 \
  --search_mode      fast \
  --scoring          vina \
  --dir              output/ \
  --verbosity        1
```

`ligands.txt` は 1 行 1 PDBQT パス。`--search_mode` は `fast`/`balance`/`detail` の 3 段階。

### 必須引数（最小セット）

| 引数 | 意味 | 必須 |
|---|---|---|
| `--receptor` | タンパク質 PDBQT | Y |
| `--ligand_index` or `--gpu_batch` | リガンドリスト or 直接指定 | Y（どちらか） |
| `--center_x/y/z` | グリッド中心 | Y |
| `--size_x/y/z` | グリッドサイズ | Y |
| `--dir` | 出力ディレクトリ | Y |
| `--scoring` | スコア関数 (`vina`/`vinardo`) | N（default=vina） |
| `--exhaustiveness` | 探索強度 | N（default=8） |

---

## 3. AutoDock Vina との差異

| 観点 | AutoDock Vina | Uni-Dock | 影響 |
|---|---|---|---|
| 実行基盤 | CPU（SIMD） | **GPU (CUDA)** | 必須: CUDA 対応 GPU（CC >= 7.0 推奨） |
| 並列単位 | 1 receptor × 1 ligand / プロセス | **1 receptor × バッチ ligand**（数百〜数千 ligand/GPU job） | Phase 2 の `dock_one_protein()` 構造と親和性高 |
| スコア関数 | Vina / Vinardo | Vina / Vinardo（AD4 系は未対応） | 互換（Phase 0-2 と同スコア出せる） |
| 速度（論文値） | 基準 | 〜**1000×**（RTX 3090 vs 単一CPUコア） | Phase 4 規模で現実的時間に収まる |
| ポーズ出力 | PDBQT | PDBQT（同形式） | 既存 `pdbqt_to_sdf` 再利用可 |
| ランダム性 | `--seed` 明示 | `--seed` 明示 | hash 決定論性は同等 |
| ADFRsuite | 不要 | 不要 | 既存 `obabel -xr` 代替で維持 |
| 追加パラメータ | なし | `--search_mode`, `--gpu_batch`, `--refine_step` | ScreeningRunner 追加引数で吸収 |

**結論**: スコア関数互換性とファイル形式互換性により、Phase 0-2 の結果と Phase 3 以降の結果を **同一 HDF5 空間で混在** できる（`source="vina"` vs `"unidock"` で区別）。

---

## 4. 既存 `DockingToolABC` との整合性検証

### 現行 ABC（`docking_automation/docking/docking.py`）

```python
class DockingToolABC(ABC):
    @abstractmethod
    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein: ...
    @abstractmethod
    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet: ...
    @abstractmethod
    def dock(self, parameters: DockingParameters) -> List[DockingResult]: ...
    @abstractmethod
    def run_docking(
        self, protein: Protein, compound_set: CompoundSet, grid_box: GridBox, ...
    ) -> DockingResultCollection: ...
```

### UniDockDocking 実装方針

| 責務 | Vina 実装 (`AutoDockVina`) | Uni-Dock 実装（新設 `UniDockDocking`） |
|---|---|---|
| `_preprocess_protein` | 水素削除 + obabel -xr → PDBQT | **同じ**（`MoleculeConverter.protein_to_pdbqt` 共有） |
| `_preprocess_compound_set` | Meeko + RDKit → PDBQT/件 | **同じ**（`compound_to_pdbqt` 共有） |
| `dock` | 1 ligand ずつ `vina.set_ligand_from_file` | **バッチ**: `ligand_index.txt` 生成 → `unidock` CLI サブプロセス呼び出し |
| `run_docking` | `dock` を多回呼び出し | `dock` を **1 回**（複数 ligand を1バッチ）で呼び出し |
| 戻り値 | `DockingResultCollection`（既存） | `DockingResultCollection`（同じ） |

**結論**: **既存 ABC シグネチャを変更せずに UniDockDocking を追加可能**。`dock` の内部実装だけがバッチ化されるが、外部呼び出し側（ScreeningRunner）にとっては同じ API。ScreeningRunner の `_dock_one_pair()` を `_dock_one_protein_batch()` に昇格させる Wave 2 変更が最小限の修正で済む。

### Phase 2 ScreeningRunner への影響

現行 `dock_one_protein(protein, compound_indices)` は既に「1 protein × N 化合物」バッチ単位。これを UniDockDocking 実装に置き換えれば、1 Dask task = 1 GPU ジョブ = バッチ投入という自然なマッピングになる。**Phase 2 の設計が Phase 3 への移行を先取りしており、追加工数は最小**。

---

## 5. インストール方法

### 5.1 公式配布

| 方式 | コマンド例 | 備考 |
|---|---|---|
| conda | `conda install -c conda-forge unidock` | 推奨。依存解決が楽 |
| pip | `pip install unidock-tools`（CLI ラッパ） | 本体バイナリは別途 |
| source build | CMake + CUDA toolkit | GPU アーキテクチャ別に最適化可 |
| Docker | `docker pull dptechnology/unidock:latest` | 開発環境推奨 |

### 5.2 CUDA 依存

- CUDA 11.x または 12.x（GPU ドライバは 525.60+ 推奨）
- GPU Compute Capability **7.0 以上**（Volta V100 / Turing T4 / Ampere A100 / Ada RTX40xx）
- VRAM 推奨: 8GB 以上（バッチサイズ依存）

### 5.3 Dockerfile への追加案

```dockerfile
# Phase 3 以降
RUN conda install -c conda-forge -y unidock \
 && python -c "import unidock_tools; print(unidock_tools.__version__)"
```

または:
```dockerfile
FROM nvidia/cuda:12.1.1-devel-ubuntu22.04 AS unidock-builder
RUN apt-get update && apt-get install -y cmake git
RUN git clone https://github.com/dptech-corp/Uni-Dock && cd Uni-Dock/unidock \
 && mkdir build && cd build && cmake .. && make -j
# 成果物を最終イメージへコピー
```

**CI 制約**: GitHub Actions の GPU ランナーは有料枠のみ。Phase 3 以降の CI は mock unidock バイナリ + 結果ファイル fixture で通す方針が妥当（cmd_002 の mock_cli.sh と同思想）。

---

## 6. Phase 4 見積もり

### 6.1 理論 speedup

論文（JCTC 2023）値より、A100 GPU 1 枚で以下の実測:
- 受容体 1 × 化合物 1M を ~40 分（`fast` モード、exhaustiveness=1）
- 1 pair あたり ~2.4 ms（バッチ並列前提）

Phase 4 目標 2×10⁸ ペアを A100 1 枚で処理した場合:
- `2e8 × 2.4e-3 s = 4.8×10⁵ s ≈ 133 時間 ≈ 5.5 日`（exhaustiveness=1）
- 4 GPU 並列: ~1.4 日
- 8 GPU 並列: ~0.7 日

本番の exhaustiveness=8 想定（Phase 4 本番品質切替）では ~10× 長く見積もる必要あり（exhaustiveness が線形コストではないため実測要）。

### 6.2 時間外挿（現行 CPU Vina との対比）

| 構成 | 2×10⁸ ペア所要時間（exh=1） |
|---|---|
| CPU Vina 1 コア（PoC 12.15s/pair） | ~77 年 |
| CPU Vina 4 並列（実測 speedup 3.17x） | ~24 年 |
| GPU Uni-Dock (A100 × 1) | ~5.5 日 |
| GPU Uni-Dock (A100 × 4) | ~1.4 日 |

**CPU 単体では Phase 4 は実質実行不可能**。GPU 採用は必須要件。

### 6.3 GPU 不可時の代替路線

GPU が用意できない場合の縮退運用オプション:

| 案 | 内容 | スケール | 評価 |
|---|---|---|---|
| A | 化合物セット縮小（10⁴ → 10³） | 2×10⁷ ペア | 2.4 年相当。非現実的 |
| B | タンパク質セット縮小（2×10⁴ → 10³） | 10⁷ ペア | 1.2 年相当。非現実的 |
| C | HPC CPU 大量投入（1000 コア） | 2×10⁸ / 12.15s ÷ 1000 ≈ 28 日 | 実行可能だが計算資源コスト大 |
| D | AFDB 全体でなく mouse essential（~2000）＋ FDA approved（~14K） | 2.8×10⁷ ペア | 3.4 年 CPU、0.8 日 GPU。GPU 採用推奨 |

**推奨**: GPU 調達が Phase 4 の前提条件。調達不可の場合は (D) で規模を下げた運用を殿へ提案。

### 6.4 Phase 4 実装スコープ案（cmd_011+ 参考）

- `docking/unidock_docking.py` 新設（`DockingToolABC` 継承）
- `DockingParameters` に `search_mode`, `gpu_batch_size` 等の UniDock 固有パラメータを追加
- `ScreeningRunner` の `backend` 引数（既に予約済）を活用、`_dock_one_protein` を UniDock 向けに差し替え
- HPCExecutor 抽象化（Dask → Slurm/PBS 対応）を並行実装
- HDF5 スキーマは Phase 2 設計どおり `source="unidock"` で識別

---

## 7. 既存コードへの変更影響

| ファイル | 変更有無 | 備考 |
|---|---|---|
| `docking/docking.py`（ABC） | **変更なし** | シグネチャ互換 |
| `docking/docking_parameters.py` | **追加のみ**（後方互換） | `UniDockSpecificParameters` 追加 |
| `docking/autodockvina_docking.py` | **変更なし** | 既存 Vina 実装は温存 |
| `docking/unidock_docking.py` | **新規** | UniDockDocking クラス |
| `converters/molecule_converter.py` | **変更なし** | PDBQT 出力共有 |
| `docking/screening_runner.py` | **軽微**（backend 切替ロジック） | Track 009-A 設計の `backend` 引数を活用 |
| `infrastructure/repositories/hdf5_docking_result_repository.py` | **変更なし** | `source` カラム既存 |
| `.devcontainer/Dockerfile` | **追加**（CUDA + Uni-Dock） | GPU 環境依存 |
| `setup.py` / `requirements.txt` | **追加**（`unidock-tools` オプション依存） | extras_require で optional 化推奨 |

**結論**: Phase 3 は **純粋な加算実装**。Vina を壊さず UniDock を並列で提供する設計が可能。

---

## 8. 推奨アクションアイテム

### Phase 3 着手前

1. **GPU 環境調達方針の確定（殿判断）**: ローカル GPU マシン / クラウド A100 / HPC 学内 GPU のいずれを前提とするか
2. **Phase 2 Track 009-F 発見の HDF5 overhead 問題解消**: Phase 3 着手と同時または直前に Sec 3.3 の protein 束ね構造へ再構築（Phase 4 規模では死活問題）
3. **mock unidock の CI 整備方針**: GPU が CI で使えない前提で、Phase 3 実装タスクに mock バイナリ作成を含める

### Phase 3 実装タスク案（cmd_011 以降の参考）

- `subtask_011_a`: UniDockDocking クラス骨格（既存ABC継承）
- `subtask_011_b`: ligand_index.txt 生成ロジック
- `subtask_011_c`: サブプロセス呼び出し・スコア/ポーズ解析
- `subtask_011_d`: ScreeningRunner `backend="unidock"` 対応
- `subtask_011_e`: mock unidock バイナリと CI 統合
- `subtask_011_f`: 小規模 GPU 実機検証（10×100）

---

## 9. 結論

- Uni-Dock は AutoDock Vina と **入出力・スコア関数互換**、既存 ABC を変更せず統合可能
- Phase 2 の ScreeningRunner 設計（1 protein × バッチ化合物 = 1 Dask task）が Uni-Dock のバッチ API と自然に対応
- GPU 調達が前提、1000×オーダーの speedup で Phase 4 を現実時間に収める
- CI は mock で対応、本番検証は GPU 実機で実施
- Phase 3 実装は純加算で Vina を温存

**Phase 3 (cmd_011) 着手可**。ただし GPU 調達方針と HDF5 構造再設計（Sec 3.3）は着手前に殿判断要請。
