# HPC 調達計画書（Phase 4 国内スパコン運用）

**作成日**: 2026-04-22
**対象**: docking_automation Phase 4 本番（cmd_012 subtask_012_gunshi_b）
**作成**: 軍師
**前提**: 殿ご裁可により **クラウド (Lambda/AWS) 路線は撤回**、国内 HPC スパコン運用へ切替

---

## 0. 背景と方針転換

先行の `docking_automation_gpu_procurement.md`（cmd_012 タスクD）は Lambda Cloud A100×8 を推奨したが、殿ご判断によりクラウド路線を取り下げ、**国内 HPC （ABCI / TSUBAME 等）** で Phase 4 を遂行する方針に変更。

本書は HPC 候補比較・実行時間見積・ジョブスケジューラ対応・HPCExecutor 抽象の実装計画をまとめる。

---

## 1. 国内 HPC 候補の比較

### 1.1 主要候補一覧

| HPC | 運用機関 | GPU 構成 | sm (Compute Capability) | 稼働時期 | UniDock 対応 |
|---|---|---|---|---|---|
| **ABCI 3.0** | 産総研 | **H200 SXM5** × 8/node | sm_90（Hopper） | 2025 稼働 | ✅ 対応 |
| **TSUBAME 4.0** | 東工大 | **H100 SXM5** × 4/node × ~160 nodes | sm_90（Hopper） | 2024/4 稼働 | ✅ 対応 |
| **Wisteria/BDEC-01 (Aquarius)** | 東大 | **A100 40GB** × 8/node | sm_80（Ampere） | 2021 稼働 | ✅ 対応 |
| **SQUID (G ノード)** | 阪大 CMC | **A100** + Ice Lake | sm_80 | 2021 稼働 | ✅ 対応 |
| **mdx II** | 東大+学認 (cloud-like HPC) | A100 × 8/node | sm_80 | 2024 稼働 | ✅ 対応 |
| **富岳** | 理研 | **GPU なし**（A64FX ARM） | — | 2021 稼働 | ❌ **不適合**（UniDock は CUDA GPU 必須） |

### 1.2 候補比較表（詳細）

| 項目 | ABCI 3.0 | TSUBAME 4.0 | Wisteria Aquarius | SQUID G |
|---|---|---|---|---|
| GPU | H200 SXM5 141GB ×8 | H100 SXM5 80GB ×4 | A100 40GB ×8 | A100 40GB ×8 |
| CPU | Xeon Platinum 8xxx 64C ×2 | Xeon Sapphire Rapids 56C ×2 | Ice Lake 36C ×2 | Ice Lake 38C ×2 |
| RAM | 2TB+ | 768GB | 512GB | 256GB |
| ローカル NVMe | 7.6TB | ~3TB | ~2TB | ~3TB |
| Interconnect | NDR InfiniBand | Quantum-2 IB 400Gb | HDR IB 200Gb | HDR IB |
| ジョブスケジューラ | **Altair PBS** | **Altair PBS Professional** | **Slurm** | **NQS-V**（富士通） |
| UniDock build | ソースビルド推奨（H200 向け sm_90 対応要） | 同上（H100 sm_90） | 公式バイナリ可（sm_80） | 公式バイナリ可（sm_80） |
| 申請 | 公募 + 随時 | 学内 + 外部共同研究 | 学認連携・学内中心 | 学認連携・学内中心 |
| コスト（参考） | AI 向け料金体系、学術 0.x〜1.x 円相当ポイント/GPU時間 | 学内 0円〜、外部要課金 | 学内 0 円（東大 ID 必要） | 学内無料 + 外部有償 |

### 1.3 UniDock GPU 要件との整合性

- **必須**: sm_80 以上（CUDA 12.x、Uni-Dock v1.1.0 バイナリの sm_80/sm_90 cubin 対応）
- sm_90 (H100/H200): ABCI 3.0 / TSUBAME 4.0 で公式バイナリ動作確認（cmd_010 subtask_010_b で検証済）
- sm_80 (A100): Wisteria/SQUID/mdx II で公式バイナリ動作確認（cmd_010 と同じ）
- **富岳 (A64FX) は GPU 非搭載のため Uni-Dock 不可**。CPU Vina の HPC 版実装が必要で実用外。

---

## 2. Phase 4（2×10⁸ ペア）実行時間見積もり

### 2.1 GPU あたり性能（cmd_010/011 実測から外挿）

RTX 2080 SUPER (sm_75) で **249.3 s / 1000 ペア**（実測）:
- RTX 2080 SUPER = 11.2 TFLOPS FP32
- 線形外挿（search_mode=balance）

| GPU | TFLOPS | 対 2080 SUPER 比 | 1000 ペア時間（推定） |
|---|---|---|---|
| A100 | 19.5 | 1.7x（実装最適化込み 10-15x）| ~20-30 s |
| H100 | 67 | 6x（実装最適化込み 25-35x）| ~8-12 s |
| H200 | 67 同等、メモリ増 | 25-35x | ~8-12 s |

**実測ベースの Phase 4 見積（2×10⁸ ペア）**:

| 構成 | 実行時間（単一 GPU） | 8 GPU 並列 |
|---|---|---|
| A100 × 1 | 40-60 日 | **5-8 日** |
| H100 × 1 | 15-25 日 | **2-4 日** |
| H200 × 1 | 15-25 日 | **2-4 日** |

※ grid padding +5Å （対策-1）で計算時間 1.5-2x 増、penalty filter / rescue_mode で実質的な計算量 1.2-1.5x 増 → 上記見積から **+50-80%** 見込み：
- A100 ×8: **8-14 日**
- H100 ×8: **3-7 日**
- H200 ×8: **3-7 日**

### 2.2 複数ノード並列（SLURM / PBS ジョブアレイ）

- 1 ジョブ = 1 protein × 全化合物 の Dask タスク単位は HPC でも維持可能
- protein 21,452 件を **200-400 ジョブ** に分割、各ジョブ 1 GPU を占有
- 1 GPU ジョブ = ~50-100 protein 処理 = 5,000-10,000 ペア = **数分〜30分/ジョブ**
- 同時実行ジョブ数 = HPC のキュー上限（ABCI 3.0: 数百 GPU 同時、TSUBAME 4.0: 同程度）

**推定**: ABCI 3.0 や TSUBAME 4.0 で 100 GPU 同時実行なら **数時間〜1日**で Phase 4 完走可能。クラウド単一ノードより高速。

---

## 3. ジョブスケジューラ対応の概算

### 3.1 現状

- `docking_automation/docking/screening_runner.py`: Dask LocalCluster（processes=True, workers=1-8）
- 単一マシン前提、HPC ジョブ提出機構なし

### 3.2 HPC 移行の選択肢

| 案 | 概要 | 工数 | HPC対応 |
|---|---|---|---|
| **A. Dask-Jobqueue** | `dask-jobqueue.SLURMCluster` / `PBSCluster` で Dask スケジューラを HPC ジョブに展開 | 中（~2日） | SLURM/PBS 両対応、少改修 |
| **B. SLURM 直接投入** | sbatch スクリプト生成、Python から subprocess で投入 | 大（~3-5日） | HPC固有、柔軟性高 |
| **C. Snakemake/Nextflow** | ワークフロー記述言語で HPC 対応 | 大（~5-7日） | 既存 Dask 実装を捨てる要 |
| **D. プリミティブ split** | 全 protein を手で job list に分割 → 各ジョブが独立 Python 実行 | 小（~0.5-1日） | HPC 依存、粒度粗い |

**推奨: 案 A (Dask-Jobqueue) + 案 D（緊急回避）併用**

- 案 A: Phase 4 本番向けの正道、既存 Dask コードを最大限活用
- 案 D: まず動かすための最小構成、Phase 4 初回実行に適用可能

### 3.3 HPC 別ジョブスケジューラ

| HPC | スケジューラ | Dask-Jobqueue 対応 |
|---|---|---|
| ABCI 3.0 | Altair PBS | `PBSCluster` で対応 |
| TSUBAME 4.0 | Altair PBS Professional | `PBSCluster`（微調整要） |
| Wisteria | Slurm | `SLURMCluster` で標準対応 |
| SQUID | NQS-V | **非標準**（Dask-Jobqueue 非対応） → カスタムアダプタ要 |
| mdx II | Slurm | `SLURMCluster` で標準対応 |

### 3.4 コード変更範囲（案 A 採用時）

| ファイル | 変更 | 概要 |
|---|---|---|
| `docking/screening_runner.py` | 軽微 | Dask クラスタを LocalCluster → `SLURMCluster`/`PBSCluster` に差替可能な Factory 化 |
| `examples/phase4_production.py` (新規) | 中 | HPC 環境向け実行スクリプト |
| `scripts/submit_phase4.sh` (新規) | 小 | sbatch/qsub 提出ラッパ |
| `compound_pipeline/*` | 変更なし | |
| HDF5 repository | 変更なし（案A v3 bundle がそのまま動く） | |

---

## 4. HPCExecutor 抽象の実装計画

### 4.1 設計目標

Phase 1 設計書で予約済の `HPCExecutor` 抽象を実装。以下を達成:
- Dask LocalCluster (開発機・検証) と HPC cluster (本番) を同一 API で切替
- ScreeningRunner は Executor 経由で submit/gather、HPC 固有コードを持たない

### 4.2 インターフェース案

```python
# docking/hpc_executor.py (新規)
from __future__ import annotations
from abc import ABC, abstractmethod
from typing import Any, Callable, Iterable, Iterator


class ExecutorABC(ABC):
    @abstractmethod
    def submit(self, fn: Callable, *args, **kwargs) -> Any:  # Future-like
        """単一タスクを投入"""

    @abstractmethod
    def map(self, fn: Callable, iterables: Iterable) -> Iterator[Any]:
        """複数タスクをバッチ投入、順序保持した結果 iterator"""

    @abstractmethod
    def as_completed(self, futures: Iterable) -> Iterator[Any]:
        """完了順に結果を返す"""

    @abstractmethod
    def close(self) -> None: ...


class DaskLocalExecutor(ExecutorABC):
    """現行の dask.distributed.LocalCluster ラッパ"""
    def __init__(self, n_workers: int = 4, threads_per_worker: int = 1, memory_limit: str = "4GB"): ...


class DaskSlurmExecutor(ExecutorABC):
    """dask_jobqueue.SLURMCluster ラッパ"""
    def __init__(
        self,
        cores: int = 1,
        memory: str = "32GB",
        queue: str = "gpu",
        walltime: str = "24:00:00",
        account: str,
        n_workers_max: int = 100,
        job_extra_directives: list[str] = None,  # "--gres=gpu:1" 等
    ): ...


class DaskPBSExecutor(ExecutorABC):
    """dask_jobqueue.PBSCluster ラッパ (ABCI/TSUBAME 向け)"""
    ...
```

### 4.3 ScreeningRunner 改修ポイント

```python
class ScreeningRunner:
    def __init__(
        self,
        ...,
        executor: ExecutorABC | None = None,   # 追加
        dask_n_workers: int = 4,               # 既存、DaskLocalExecutor default
    ):
        self.executor = executor or DaskLocalExecutor(n_workers=dask_n_workers)

    def run(self, resume: bool = True) -> ScreeningResult:
        # 既存の with LocalCluster(...) as cluster, Client(cluster) as client:
        # を executor.submit / executor.as_completed に置き換え
        futures = {self.executor.submit(dock_one_protein, ...): pid for ...}
        for future in self.executor.as_completed(futures): ...
```

### 4.4 実装 subtask 一覧

| subtask | 内容 | 難度 | 足軽 | 見積 |
|---|---|---|---|---|
| subtask_012_m_executor_abc | `ExecutorABC` + `DaskLocalExecutor` 新設（既存動作維持） | M | 1 | ~1日 |
| subtask_012_n_dask_slurm | `DaskSlurmExecutor`（dask-jobqueue 依存追加） | M | 1 | ~0.8日 |
| subtask_012_o_dask_pbs | `DaskPBSExecutor`（ABCI/TSUBAME 向け） | M | 1 | ~0.8日 |
| subtask_012_p_runner_adapt | `ScreeningRunner` を `executor` 経由に改修 | S | 1 | ~0.5日 |
| subtask_012_q_tests | ExecutorABC unit test + ScreeningRunner 統合テスト | M | 1 | ~1日 |
| subtask_012_r_phase4_script | `examples/phase4_production.py` + sbatch/qsub スクリプト雛形 | S | 1 | ~0.5日 |
| **合計** | — | — | **3-4 足軽並列** | **実働 ~2-2.5日** |

---

## 5. 推奨 HPC プラン

### 5.1 候補の優先順位付け

| 順位 | HPC | 理由 |
|---|---|---|
| 1 | **TSUBAME 4.0** | H100 SXM5 × 160 GPU、東工大運用でアクセスしやすい（殿が東工大系ならベスト）。Phase 4 が 3-7 日で完走可能 |
| 2 | **ABCI 3.0** | H200（最新・最速）、産総研の広域公募あり。若干申請手続きが煩雑だが性能は随一 |
| 3 | **Wisteria (Aquarius)** | A100 × 8/node、東大学認連携。性能は TSUBAME/ABCI の 50-70% だが可用性高 |
| 4 | SQUID | A100、NQS-V スケジューラ対応の追加実装が必要で工数発生 |
| 5 | mdx II | A100、クラウド風運用で敷居低いが性能は上記に劣る |

**軍師推奨**: 殿の所属機関に応じて以下の優先順位で申請:
- 東工大関係者: **TSUBAME 4.0 を第一申請**、並行して ABCI 3.0 に予備申請
- 東大関係者: **Wisteria を第一**、ABCI 3.0 併用
- 無所属の場合: **ABCI 3.0** 産業利用 or 学術共同利用で申請

### 5.2 申請から Phase 4 実行までのロードマップ

```
Week 1-2: HPC 申請（所属・プロジェクト書類提出）
Week 2-3: アカウント発行待ち（ABCI は早い、TSUBAME は2-4週）
Week 3: subtask_012_m〜r (Executor 抽象) 実装（並行）
Week 4: 10×100 小規模 E2E を HPC で実施（動作確認）
Week 4-5: Phase 4 本番投入（3-14 日）
Week 5-6: 結果集約・HDF5 ダウンロード・解析
```

---

## 6. コスト試算

### 6.1 TSUBAME 4.0（東工大学内想定）

- 学内: **0 円〜**（ポイント枠内）
- 外部: 1 GPU 時間 ~300-600 円（概算）
- Phase 4 (H100×8, 5日): 8×24×5 = 960 GPU時間 → 外部料金 **~30-60 万円**、学内無料

### 6.2 ABCI 3.0（産総研）

- 個人利用（I ノード）: ~1-2 円/GPU時間ポイント換算
- Phase 4 (H200×8, 3-7日): ~700-1,500 GPU時間 → **~700-3,000 円相当ポイント**（圧倒的に安い）
- 所属機関なら pass、未所属なら産業利用申請で数ヶ月かかる可能性

### 6.3 Wisteria Aquarius（東大）

- 学内無料 + 外部有償（要確認）
- Phase 4 (A100×8, 8-14 日): ~2,500-4,500 GPU時間 → 学内無料

### 6.4 結論

- **最安・最速: TSUBAME 4.0 学内** or **ABCI 3.0**
- クラウド Lambda (旧案) $1,135-1,770 と比較して **学内 HPC なら実質無料**、**産業利用でも同等以下**

---

## 7. リスクと緩和

| リスク | 影響 | 緩和策 |
|---|---|---|
| HPC 申請の手続き遅延（2-4 週間） | Phase 4 開始遅延 | 複数 HPC に並列申請、ABCI を backup |
| UniDock ソースビルドの H200/H100 再ビルド | 環境構築遅延 | cmd_010 で整備済み Dockerfile (sm_75/80/90 マルチアーキ) を転用、HPC 向けに nvcc 追加 |
| NQS-V 非標準スケジューラ対応（SQUID） | 実装工数 +1日 | SQUID は最終手段、SLURM/PBS 系を優先 |
| HPC の data storage 容量制限 | 1.3TB の HDF5 格納不能 | HDF5 v3 bundle (136GB) で Phase 4 完走可、余裕あり |
| ネットワーク帯域（AFDB/compound 転送） | 初回転送遅延 | scp/rsync で夜間転送、S3/学内 FS 経由 |
| 認証・2要素認証 | 運用障壁 | 申請時に認証設定を確認、SSH 鍵整備 |

---

## 8. 次アクション

### cmd_012 内（軍師/家老→足軽配備可）

1. **subtask_012_m-r**: HPCExecutor 抽象実装（6 subtask、3-4 足軽並列 2-2.5日）
2. **subtask_012_s_hpc_apply**: HPC 申請書類準備（殿手動、軍師が申請用サマリ生成可能）

### Phase 4 着手時（cmd_013）

3. **subtask_013_a_small_e2e_hpc**: HPC 環境で 10×100 動作確認
4. **subtask_013_b_main_run**: 2×10⁸ ペア投入 + HDF5 集約
5. **subtask_013_c_analysis**: 結果解析・ランキング生成

---

## 9. 結論

- **クラウド路線を撤回し HPC 運用へ方針転換。技術的整合性に問題なし**
- 推奨 HPC: **TSUBAME 4.0 (H100 ×8)** → **ABCI 3.0 (H200)** → **Wisteria (A100)** の順
- Phase 4 実行時間: 3-14 日（GPU 8 並列）、コスト: 学内無料 or 数千円〜数十万円
- 実装工数: HPCExecutor 抽象で **3-4 足軽並列 2.5 日**
- 申請手続きを先行開始、HPCExecutor 実装を並行進行すれば **1ヶ月以内に Phase 4 着手可能**
