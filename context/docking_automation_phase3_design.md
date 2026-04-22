# Phase 3 設計書: UniDockDocking + ScreeningRunner バックエンド切替

**作成日**: 2026-04-22
**対象**: docking_automation v2 Phase 3（cmd_011）
**作成**: 軍師（subtask_011_design_estimate）
**前提**: cmd_010 完全完了（Uni-Dock v1.1.0 導入済 / GPU E2E PASS 5/5 差 0.018 kcal/mol / libcuda fix 適用済 / cuda-compat 削除済 commit c4fe84c）

---

## 1. Phase 3 スコープ

| 含む | 含まない |
|---|---|
| `UniDockDocking(DockingToolABC)` 新規実装 | HPC Executor 抽象化（Phase 4） |
| バッチ ligand 投入（--ligand_index） | マルチ GPU 並列（Phase 4） |
| ScreeningRunner への backend 切替 | HDF5 構造再設計（Phase 4）|
| 10×100 GPU E2E 検証（開発機 RTX 2080 SUPER） | 本番 2×10⁸ 実行（Phase 4） |
| A100/H100 への移植確認計画 | A100/H100 実機ベンチ |

Phase 3 完了基準: **ScreeningRunner(backend="unidock") で 10×100 ペアを GPU 実行**、Vina 結果との score 相関 > 0.9、resume 冪等性維持。

---

## 2. UniDockDocking クラス設計

### 2.1 責務

```
UniDockDocking (DockingToolABC 継承)
 ├─ _preprocess_protein(protein)     # AutoDockVina と共通実装 (obabel -xr)
 ├─ _preprocess_compound_set(cset)   # AutoDockVina と共通実装 (Meeko)
 ├─ dock(parameters)                 # CLI subprocess 呼出し・バッチモード
 └─ run_docking(protein, cset, gb, ...)  # dock を 1 回呼んで DockingResultCollection 返却
```

### 2.2 I/F（型ヒント付き疑似コード）

```python
# docking_automation/docking/unidock_docking.py
from __future__ import annotations
import subprocess, tempfile, shutil, json
from pathlib import Path
from typing import List, Optional
from docking_automation.docking.docking import DockingToolABC
from docking_automation.docking.docking_parameters import (
    CommonDockingParameters, DockingParameters, SpecificDockingParametersABC,
)
from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet


class UniDockParameters(SpecificDockingParametersABC):
    """Uni-Dock 固有パラメータ (SpecificDockingParametersABC の実装)。"""
    def __init__(
        self,
        search_mode: str = "fast",       # fast | balance | detail
        scoring: str = "vina",           # vina | vinardo
        num_modes: int = 3,
        max_step: Optional[int] = None,  # override for detail
        refine_step: Optional[int] = None,
        seed: int = 1,                   # 決定論性確保 (Phase 3 で追加)
        verbosity: int = 1,
    ) -> None: ...


class UniDockDocking(DockingToolABC):
    """Uni-Dock v1.1.0 をサブプロセス呼び出しで利用する DockingTool 実装。

    AutoDockVina との主な差異:
    - 1 run = 1 protein × バッチ ligand（--ligand_index で列挙）
    - GPU 実行、1000× 級 speedup
    - Meeko + obabel 前処理は共通
    """
    UNIDOCK_BINARY: str = "unidock"  # PATH 前提、compat でフルパス化可能

    def __init__(self, binary_path: Optional[str] = None) -> None: ...

    # --- 前処理: AutoDockVina と同一実装を継承 or コピー ---
    def _preprocess_protein(self, protein):
        # 既存の MoleculeConverter.protein_to_pdbqt を呼ぶだけ
        ...

    def _preprocess_compound_set(self, compound_set):
        # 既存の MoleculeConverter.compound_to_pdbqt を呼ぶだけ
        ...

    # --- 中核: dock (バッチ) ---
    def dock(self, parameters: DockingParameters) -> List[DockingResult]:
        """1 protein × N 化合物を 1 回の unidock CLI 呼び出しで処理する。
        戻り値は N 個の DockingResult。
        """
        self._validate_params(parameters)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp = Path(tmpdir)
            lig_index = tmp / "ligands.txt"
            out_dir = tmp / "out"
            out_dir.mkdir()

            # 1) ligand_index.txt 生成
            ligand_paths = parameters.preprocessed_compounds.file_paths
            lig_index.write_text("\n".join(str(p) for p in ligand_paths) + "\n")

            # 2) unidock CLI 実行
            cmd = self._build_cmd(parameters, lig_index, out_dir)
            completed = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if completed.returncode != 0:
                # 失敗時は全ペアを error=... で返す
                return self._make_error_results(parameters, completed.stderr)

            # 3) ポーズ・スコア解析
            results = self._parse_outputs(parameters, out_dir, completed.stdout)
            return results

    # --- run_docking: Phase 2 既存 API と一致 ---
    def run_docking(
        self, protein, compound_set, grid_box,
        parameters=None, repository=None, **kwargs,
    ) -> DockingResultCollection:
        prep_p = self._preprocess_protein(protein)
        prep_c = self._preprocess_compound_set(compound_set)
        params = self._build_docking_parameters(prep_p, prep_c, grid_box, parameters)
        results = self.dock(params)
        # repository への保存は ScreeningRunner 側責務（Phase 2 と同じ）
        return DockingResultCollection.from_list(results)

    # --- 内部実装 ---
    def _build_cmd(self, params, lig_index, out_dir) -> List[str]:
        p = params.specific_parameters   # UniDockParameters
        return [
            self.UNIDOCK_BINARY,
            "--receptor", str(params.preprocessed_protein.file_path),
            "--ligand_index", str(lig_index),
            "--center_x", str(params.grid_box.center[0]),
            "--center_y", str(params.grid_box.center[1]),
            "--center_z", str(params.grid_box.center[2]),
            "--size_x", str(params.grid_box.size[0]),
            "--size_y", str(params.grid_box.size[1]),
            "--size_z", str(params.grid_box.size[2]),
            "--scoring", p.scoring,
            "--search_mode", p.search_mode,
            "--num_modes", str(p.num_modes),
            "--seed", str(p.seed),
            "--dir", str(out_dir),
            "--verbosity", str(p.verbosity),
        ]

    def _parse_outputs(self, params, out_dir, stdout) -> List[DockingResult]:
        """out_dir/{ligand_stem}_out.pdbqt と stdout からスコアを抽出。"""
        # Uni-Dock は stdout に "REMARK VINA RESULT:" 形式でスコアを出す
        # 各 ligand について 1 つの top-pose を取り出して DockingResult 構築
        ...
```

### 2.3 DockingParameters への影響

- 既存 `CommonDockingParameters` に `exhaustiveness` があるが、Uni-Dock は `exhaustiveness` を直接取らず `search_mode` (fast/balance/detail) で制御する → UniDockParameters 側で `search_mode` を持つ
- 既存 `AutoDockVinaParameters` は変更なし
- **共通化提案**: `CommonDockingParameters.exhaustiveness` を `UniDockDocking._to_search_mode(exhaustiveness)` で fast/balance/detail にマッピングし、呼び出し側 API 統一（ScreeningRunner が `exhaustiveness` だけ指定すれば両 backend 動作）

### 2.4 決定論性確保

- `--seed N` を明示指定（UniDockParameters.seed=1 デフォルト）
- Uni-Dock CUDA 実装は atomic ops 順序依存の非決定性がわずかにあるが、seed 固定で通常 ±0.01 kcal/mol 程度の揺らぎに収まる（Phase 2 で測定した差 0.018 kcal/mol と整合）
- content_hash は pose_blob (gzip SDF) ベースのため、座標 ±0.001Å のブレで hash 不一致の可能性 → Phase 3 E2E で hash drift 発生率を実測し、許容範囲を殿判断

---

## 3. ScreeningRunner バックエンド切替機構

### 3.1 現状確認（Phase 2）

Phase 2 実装（`screening_runner.py`）では `dock_one_protein()` が module-level 関数で、内部で `AutoDockVina()` を直接生成している（cmd_009 Track C の cloudpickle 対応のため）。これを backend 切替可能にする。

### 3.2 改修方針

**案 A（推奨）: DockingTool ファクトリ経由**

```python
# screening_runner.py
class ScreeningRunner:
    def __init__(
        self,
        ...,
        backend: str = "vina",   # "vina" | "unidock"
        backend_params: Optional[dict] = None,
    ) -> None:
        self.backend = backend
        self.backend_params = backend_params or {}

    def _make_docking_tool(self):
        if self.backend == "vina":
            from docking_automation.docking.autodockvina_docking import AutoDockVina
            return AutoDockVina(**self.backend_params)
        elif self.backend == "unidock":
            from docking_automation.docking.unidock_docking import UniDockDocking
            return UniDockDocking(**self.backend_params)
        raise ValueError(f"unknown backend: {self.backend}")
```

module-level `dock_one_protein()` も backend 引数を追加:
```python
def dock_one_protein(protein_info, compound_info, grid_box, backend: str, backend_params: dict):
    tool = _make_docking_tool(backend, backend_params)
    return tool.run_docking(...)
```

**利点**:
- 既存 `DockingToolABC` を変更せず、呼び出し側だけで切替可能
- HDF5 には `source="vina"` / `"unidock"` と記録（Phase 2 スキーマ既対応）
- 単体テストは backend モックで両系統を網羅可能

### 3.3 後方互換

- `backend` 引数なしで呼ぶと `"vina"` 既定 → Phase 2 E2E は完全無影響
- 既存テスト（test_screening_runner.py 7件）は `backend="vina"` 暗黙で動くはず、必要なら明示的に指定

---

## 4. GPU E2E テスト設計

### 4.1 土台

`examples/unidock_e2e_test.py`（足軽2号 commit 0426635 で新設）が最小スモーク。Phase 3 ではこれを **ScreeningRunner 経由** に格上げ。

### 4.2 新テストスクリプト案

`examples/phase3_gpu_e2e.py`:
```python
from docking_automation.docking.screening_runner import ScreeningRunner
from docking_automation.molecule.protein_set import ProteinSet
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.docking.grid_box_cache import GridBoxCache

# 10 protein × 100 compound = 1000 ペア
protein_set = ProteinSet.from_directory("examples/input/afdb_mouse")  # 10 件
compound_set = CompoundSet("examples/input/ALDR/actives_final.sdf.gz")[:100]
grid_box_cache = GridBoxCache.from_file("cache/grid_boxes.json")

runner = ScreeningRunner(
    protein_set, compound_set, grid_box_cache,
    hdf5_path="examples/output/phase3_unidock_hdf5/docking_results.hdf5",
    backend="unidock",
    backend_params={},
    dask_n_workers=1,  # GPU は単一、Dask は進捗ログ用に形式的に
)
result = runner.run(resume=True)
```

### 4.3 検証項目

1. **動作性**: 1000 ペア完了、failed=0、elapsed_sec 記録
2. **性能**: Vina (Phase 2 E2E) 比で >10× speedup（RTX 2080 SUPER 期待値）
3. **スコア相関**: Vina 同ペアと Pearson r > 0.9（個別 ±0.02 kcal/mol 程度を許容）
4. **冪等性**: 2回目実行で reused=1000 / new=0
5. **HDF5 content_hash**: 1回目と2回目で同一

### 4.4 追加 unit test

- `tests/docking/test_unidock_docking.py`: UniDockDocking 単体（mock CLI）
- `tests/docking/test_screening_runner.py` 拡張: `backend="unidock"` でも既存テストが通る（mock dock）

---

## 5. Phase 3 実装見積もり

### 5.1 subtask 分割

| subtask | 内容 | 難度 | 足軽 | 見積 |
|---|---|---|---|---|
| subtask_011_a | UniDockDocking 骨格（`dock_one_protein` module-level + backend factory） | M | 1 | ~1日 |
| subtask_011_b | UniDockParameters + CLI builder + output parser | M | 1 | ~1日 |
| subtask_011_c | ScreeningRunner backend 切替（案 A） | S | 1 | ~0.5日 |
| subtask_011_d | unit test: test_unidock_docking.py（mock CLI） | M | 1 | ~0.7日 |
| subtask_011_e | unit test: ScreeningRunner backend 切替テスト | S | 1 | ~0.3日 |
| subtask_011_f | examples/phase3_gpu_e2e.py + 10×100 実機 E2E | M | 1 | ~0.7日（+ GPU時間） |
| subtask_011_g | Phase 3 実装報告書 + スコア相関 Vina vs UniDock 結果 | S | 1 | ~0.3日 |
| **合計** | — | — | **3足軽並列** | **実働 1.5-2日** |

### 5.2 並列計画

- Wave 1（並列3）: subtask_011_a / 011_b / 011_c
- Wave 2（並列2）: subtask_011_d / 011_e（Wave 1 依存）
- Wave 3（単一）: subtask_011_f（統合）→ 011_g

### 5.3 Phase 2 からの差分

| 変更対象 | 影響 |
|---|---|
| `docking/unidock_docking.py` | **新規** |
| `docking/docking_parameters.py` | **追加のみ** (UniDockParameters) |
| `docking/screening_runner.py` | **軽微** (backend 引数と factory 追加) |
| `docking/docking.py` (ABC) | **変更なし** |
| `docking/autodockvina_docking.py` | **変更なし** |
| HDF5 スキーマ | **変更なし** (`source` カラム既対応) |
| `compound_pipeline/` | **変更なし** |
| `ProteinSet` / `GridBoxCache` | **変更なし** |

**結論**: Phase 3 は純加算実装で、Phase 2 テスト（157 passed）を破らずに追加可能。

---

## 6. Phase 4 本番への考慮

### 6.1 A100/H100 へのポータビリティ

**GPU architecture**:
- A100: sm_80（Ampere）
- H100: sm_90（Hopper）
- 開発機 RTX 2080 SUPER: sm_75（Turing）

**現状の Uni-Dock バイナリ（sm_75 patch 済）での対応**:
- 前提: cmd_010 で `CMAKE_CUDA_ARCHITECTURES="75;80;90"` マルチアーキビルドが適用されているか確認（未確認なら次 subtask）
- 適用済みなら A100/H100 でも同一 Docker イメージで動作

**推奨アクション**:
- Phase 3 subtask_011_a 開始前に Dockerfile の CMAKE_CUDA_ARCHITECTURES 設定を確認
- sm_75 のみビルドだった場合、`subtask_011_x_multi_arch_build` を別途起票

### 6.2 クラウド GPU 環境（殿ご裁可 Q1-A）

- **候補**: AWS EC2 p4d (A100) / Google Cloud A2 / Azure ND A100 v4
- **デプロイ**: 本プロジェクト Docker イメージをクラウド GPU インスタンスで起動
- **実行**: ScreeningRunner を `backend="unidock"` で大規模実行
- **データ転送**: AFDB pdb (6.7GB) + compound SDF を事前アップロード、HDF5 結果をローカル回収

**Phase 4 subtask 起票候補（Phase 3 完了後）**:
- `subtask_012_a`: クラウド GPU 環境 Terraform / Ansible 化
- `subtask_012_b`: HDF5 sharding（protein 束ね構造）で overhead 解消
- `subtask_012_c`: HPCExecutor 抽象化（Dask → SLURM）
- `subtask_012_d`: 本番 2×10⁸ 投入スクリプト
- `subtask_012_e`: 結果集約・解析レポート

### 6.3 性能見積もりの再確認

Phase 3 RTX 2080 SUPER 実測（subtask_011_f 完了後に確定）:
- 10×100 = 1000 ペアを exh=fast で ~数分想定
- 1000 ペア / 時間 → A100 で 10-20× 速く → 2×10⁸ ペアが **A100×1 で 3-5 日**見込み

軍師 Uni-Dock 調査レポート Sec 6.2 の「A100×1 で ~5.5 日」は exh=1 CPU Vina 換算値から理論外挿。実測ベースで Phase 4 計画を更新。

---

## 7. Open Questions（殿判断事項）

| # | テーマ | 軍師推奨 |
|---|---|---|
| Q1 | backend パラメータを ScreeningRunner `__init__` に追加（案 A）で良いか | 採用 |
| Q2 | UniDockParameters.seed デフォルト値（固定 1 or ランダム生成） | **固定 1**（決定論性優先） |
| Q3 | score 相関しきい値（Vina vs UniDock） | **Pearson r > 0.9** |
| Q4 | sm_75+sm_80+sm_90 マルチアーキビルドを Phase 3 着手前に先行実施するか | **先行実施推奨**（本番移行時の手戻り防止） |
| Q5 | Phase 3 E2E を Dask LocalCluster 経由にするか単一プロセスで直接呼ぶか | **Dask 経由**（Phase 4 への接続性保持、ただし workers=1） |

---

## 8. 完了基準（Phase 3）

1. `UniDockDocking(DockingToolABC)` + `UniDockParameters` 実装
2. `ScreeningRunner(backend="unidock")` が動作し既存 test 回帰なし
3. `examples/phase3_gpu_e2e.py` が 10×100 = 1000 ペア完了、failed=0
4. Vina vs UniDock score 相関 Pearson r > 0.9
5. 冪等性（2回目で reused=1000 / new=0）維持
6. pytest SKIP=0（Phase 2 範囲内、Phase 3 新規テスト含む）
7. Phase 3 実装報告書 commit + push

---

**設計確定待ち**: Q1-Q5 の殿ご裁可を受けて本設計書を最終化、cmd_011 subtask 群を起票。
