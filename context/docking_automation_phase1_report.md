# Phase 1 実装報告書

**作成日**: 2026-04-22
**対象**: docking_automation v2 ブランチ Phase 1（cmd_007）
**作成**: 軍師（subtask_010_gunshi Part B）

---

## 1. 実装サマリ

Phase 1 は 4 トラック並列で実装。全トラック commit + push 完了、単体テスト全件 PASS、pytest 全体 SKIP=0（Phase 2 範囲内で）。

| トラック | 内容 | 担当 | commit | 実装規模 | テスト |
|---|---|---|---|---|---|
| A | ProteinSet | 足軽1号 | `b5fd8c0` | protein_set.py 123行 + test 168行 | 19 passed / 0 skipped |
| B | GridBoxCache | 足軽2号 | `5038c05` | grid_box_cache.py 158行 + test 156行 | 12 passed / 0 skipped |
| B2 | GridBoxCache.get_with_policy 追加（Phase 2 cmd_009 で拡張） | 足軽4号 | `68bb6c2` | +31行 + test +46行 | 17 passed / 0 skipped |
| C | AFDB マウス v6 bulk DL スクリプト | 足軽3号 | `9318932` | scripts/preprocess_afdb_mouse.py 165行（6関数） | DL 成功 (3.6GB / 2026-04-22 11:42:54) |
| D | 化合物前処理パイプライン | 足軽4号 | `0544557` | compound_pipeline/preprocess.py 112行 + test 66行 + golden 3件 | 6 passed / 0 skipped |

---

## 2. トラック別詳細

### Track A — ProteinSet（集約ルート）

- 実装ファイル: `docking_automation/molecule/protein_set.py` (123 lines)
- テスト: `tests/molecule/test_protein_set.py` (168 lines)
- **I/F レビュー指示への準拠**:
  - `from __future__ import annotations`
  - `Protein(path)` コンストラクタ直接呼び出し（`Protein.create()` の UUID 生成を回避）
  - `content_hashes()` キャッシュ実装（`cached` dict）
  - `from_directory(root)` / `from_afdb_mouse(root, max_residues=2000, limit=None)` ファクトリ
  - 重複 `protein_id` の `__init__` 検出による fail-fast
  - AF-*.pdb glob + 残基数 > 2000 フィルタ

### Track B — GridBoxCache（JSON 永続キャッシュ）

- 実装ファイル: `docking_automation/docking/grid_box_cache.py` (158 lines → 189 after cmd_009 Track 009-D)
- テスト: `tests/docking/test_grid_box_cache.py` (156 → 202 lines)
- **I/F レビュー指示への準拠**:
  - `TYPE_CHECKING` guard で `ProteinSet` 遅延参照 → 循環 import 完全回避
  - `GridBoxCacheEntry.to_dict/from_dict` で np.ndarray ↔ list[float] 変換を一元化
  - `missing_ids(protein_set)` で `protein_content_hash` 不一致も検出（AFDB 更新時の自動無効化）
  - `atomic_save`（tempfile + os.replace）
- **Phase 2 で追加**: `get_with_policy(protein_id, missing_policy, fallback_center, fallback_size)`
  - `skip` / `error` / `fallback_centroid` の 3 ポリシーを実装
  - 既存 `get()` は後方互換維持

### Track C — AFDB マウス v6 bulk 前処理スクリプト

- 実装ファイル: `scripts/preprocess_afdb_mouse.py` (165 lines)
- 6 関数構成:
  1. `preflight(required_gb)` — ディスク容量検証（21TB 空き / 要 15GB OK）
  2. `download(resume=True)` — `wget --continue --tries=3 --timeout=60`
  3. `validate_tar()` — tar 整合性チェック
  4. `extract(skip_existing=True)` — `*.pdb.gz` のみ `extracted/` へ展開
  5. `decompress_pdb(limit=None)` — `pdb.gz → pdb`、limit で件数制限可
  6. `validate(sample_n=10)` — PDB サンプルパース検証
- `.gitignore` に `data/afdb/` を追記
- **DL 完了**: `UP000000589_10090_MOUSE_v6.tar` 3.6 GB / 2026-04-22 11:42:54 取得成功
- **命名規則**: `AF-{uniprot}-F1-model_v6.pdb`（設計書・ProteinSet.from_afdb_mouse() と完全整合）

### Track D — 化合物前処理パイプライン

- 実装ファイル: `docking_automation/compound_pipeline/__init__.py` + `preprocess.py` (112 lines)
- テスト: `tests/pipeline/__init__.py` + `test_compound_hash_stability.py` (66 lines)
- golden hash セット: `tests/pipeline/golden_hashes.json` (3 件: aspirin / ethanol / ibuprofen)
- **仕様準拠**（context/docking_automation_compound_pipeline.md）:
  1. Dimorphite-DL 2.0.2 `protonate_smiles(..., max_variants=1)` — ルールベース決定論的
  2. RDKit `MolStandardize.TautomerEnumerator.Canonicalize()`
  3. ETKDGv3 + `params.randomSeed = 42`
  4. `AllChem.MMFFOptimizeMolecule(maxIters=2000)`
  5. `Chem.AddHs(mol, addCoords=True)`
  6. `SDWriter` で固定精度出力 → content_hash 対象
- 決定論テスト (`test_same_smiles_same_sdf`) + クロスセッション golden hash テスト（3件）全 PASS

---

## 3. テスト結果サマリ

| テストスイート | passed | skipped | failed | error |
|---|---|---|---|---|
| test_protein_set.py | 19 | 0 | 0 | 0 |
| test_grid_box_cache.py | 17 | 0 | 0 | 0 |
| test_compound_hash_stability.py | 6 | 0 | 0 | 0 |
| Phase 1 合計 | **42** | **0** | **0** | **0** |

pytest 全体（Phase 2 完了時点）: **157 passed / 31 skipped / 0 failed / 0 error**
- 残留 SKIP 31 件は Phase 3/4 スコープの未実装プレースホルダー（殿ご裁可 A により Phase 2 範囲での SKIP=0 と定義）

---

## 4. Phase 2 連携確認

Phase 2 実装（cmd_009）にて Phase 1 成果物の統合動作が確認された:

- **ProteinSet**: ScreeningRunner のイテレーション入力として採用、`content_hashes()` を差分検出に使用
- **GridBoxCache**: `get_with_policy(policy="skip")` で未登録 protein を除外、事前 fpocket 計算不要
- **compound_pipeline**: `preprocess.py` で preprocessed SDF 生成、ScreeningRunner 入力として content_hash 安定化
- **E2E 10×10 検証**（commit `224412b`, 足軽2号 subtask_009_f）:
  - fpocket で全 10 protein の GridBox 自動計算成功
  - Run1: 100 ペア全件ドッキング成功（failed=0, 162s）
  - Run2 (resume=True): reused=100 / new=0 / 0.06s → **冪等性完全検証 OK**

---

## 5. 残課題・注記

### Track C — tar 解凍・pdb 整備（進行中）

本報告書作成時点（2026-04-22 12:40）で `data/afdb/pdb/` は空。`data/afdb/raw/UP000000589_10090_MOUSE_v6.tar` (3.6 GB) は取得済。足軽3号が並行して `extract → decompress_pdb → validate` を実施中。完了後に pdb ファイル数・validate 結果を補足報告予定。Phase 2 E2E は既存 `examples/input/afdb_mouse/` の 10 protein サンプルで動作検証済のため、AFDB 全体展開は Phase 3 本番実行時の前提として完了を要する。

### Phase 4 サイズリスク（cmd_009 Track 009-F で発見）

E2E 10×10 実測で HDF5 **bytes_per_pair = 7189.5B**（Phase 2 設計見積 ~850B の約 8 倍）。Phase 4 外挿すると **1339 GB** となり目標 200 GB を大きく超過。原因は HDF5 グループ/データセット overhead (~6.5KB/group)。Phase 2 設計書 Sec 3.3 にて示した「protein ごとに 1 データセット、化合物を行として束ねる」構造への再構築を Phase 3/4 着手前に強く推奨。

### 副次バグ修正

- Track 009-C の `dock_one_protein` が `MoleculeConverter.sdf_to_pdbqt` 未実装メソッドを呼んでおり、Track 009-F (足軽2号) が `_convert_to_pdbqt_with_meeko` をラップするメソッドを追加して解消（実機 Vina 動作が初めて成立）。

---

## 6. 結論

Phase 1 の 4 トラック全て実装完了。設計書（`context/docking_automation_phase1_design.md`）への準拠、I/F レビュー指示への反映、Phase 2 連携確認（E2E）まで完遂。Phase 3（Uni-Dock 統合）および Phase 4（HPC Executor 抽象 + HDF5 構造再設計）へ進める状態。

**設計確定日**: 2026-04-22 / 殿ご裁可による Phase 2 完了基準緩和（SKIP=0 は Phase 2 範囲内）を反映
