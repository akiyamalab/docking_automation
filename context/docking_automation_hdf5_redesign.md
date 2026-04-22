# HDF5 スキーマ再設計案（Phase 4 overhead 解消）

**作成日**: 2026-04-22
**対象**: docking_automation Phase 4（cmd_012 subtask_012_gunshi_a タスクC）
**作成**: 軍師

---

## 1. 現行スキーマの問題

### 実測値（Phase 2 Track 009-F / Phase 3 011_k）

| 指標 | 実測 |
|---|---|
| Phase 2 E2E (10×10=100 ペア) | 702 KB |
| Phase 3 E2E (10×100=1000 ペア, 有効690) | 6.7 MB (拡張後) |
| bytes/pair (Phase 2) | 7,189 B |
| bytes/pair (Phase 3 corrected) | ~7,000 B |

### Phase 4 外挿

```
2×10⁸ ペア × 7,189 B = 1,437,800 MB = 1,338 GB
```

**目標 200 GB を 6.7 倍超過**。主因は HDF5 のグループ/データセット overhead（~6,500 B/group）。

### 構造的原因

現行:
```
/results/
  {protein_hash_64char}/              ← グループ overhead ~400B
    {compound_hash_64char}/           ← グループ overhead ~400B
      docking_score        (f32)      ← scalar dataset overhead ~1KB
      pose_blob            (bytes)    ← dataset overhead + 圧縮pose 680B
      computed_at          (str)      ← scalar dataset overhead ~500B
      source               (str)      ← scalar dataset overhead ~500B
      top_n                (i32)      ← scalar dataset overhead ~500B
```

1 pair あたり 5 scalar dataset + 2 group = ~6.5 KB の固定 overhead。実データは ~680 B のみ。**圧縮率 10%**。

---

## 2. 再設計案の比較

### 案 A: Protein 束ね構造（chunked 2D arrays）

```
/results/
  {protein_hash}/                        ← 1 protein = 1 group
    compound_hashes   (S64, shape=(N,))  ← 文字列配列
    docking_scores    (f32, shape=(N,))  ← スコア1D配列
    pose_blobs        (vlen bytes, shape=(N,))  ← 可変長bytes
    computed_at       (S24, shape=(N,))  ← ISO8601タイムスタンプ
    source            (S8,  shape=(N,))  ← "vina"/"unidock"
    top_n             (i32, shape=(N,))  ← default 1
```

- 1 protein 毎に **1 group + 5 datasets**（overhead ~6.5KB）で **全化合物 N 件**を保持
- overhead は protein 数のみに依存: 2×10⁴ × 6.5 KB = **130 MB** (現行比 1万倍削減)
- 実データ: 2×10⁸ × 680 B = 136 GB

**合計**: ~136.1 GB ≤ 目標 200 GB ✅（62% 余裕）

**アクセスパターン**:
- 取得: `hdf5['/results/{protein_hash}']['compound_hashes'][:]` で indices 探索 → `['docking_scores'][idx]`
- 差分検出: 起動時に全 compound_hashes を set 化して cache（メモリ ~6GB 想定）
- 書込: worker が 1 protein ぶんの結果を返す → メインで appendモード or 上書き

**評価**:
| 観点 | 評価 |
|---|---|
| サイズ | ◎ 136 GB ≤ 200 GB |
| 実装容易さ | ○ hdf5_docking_result_repository.py を中〜大改修 |
| アクセス速度 | ○ 1 protein 単位 bulk read が速い |
| Resume 差分検出 | ○ compound_hashes 全ロード後に set 比較 |
| 並列書込み | △ 既存の「メイン専有 write」パターンを維持、ただし worker 返却単位を protein 単位に揃える必要 |

### 案 B: UniProt 先頭2文字シャーディング

```
data/hdf5/
  AF-A0.h5   ← AF-A0* 全 protein
  AF-A1.h5
  ...
  AF-Q9.h5
  (~300 shards 想定)
```

- 各シャード内は案 A と同じ構造
- 複数プロセスから独立シャードに並列書込可能（Phase 4 でのスループット向上）
- 単一ファイル肥大化リスクを回避

**評価**:
| 観点 | 評価 |
|---|---|
| サイズ | ◎ 案Aと同じ |
| 実装容易さ | △ shard routing ロジック追加 |
| アクセス速度 | ○ shard単位で並列読み可能 |
| Resume 差分検出 | △ shard 単位の set ロードで分散 |
| 並列書込み | ◎ shard ごとに独立 writer 可能 |

### 案 C: SQLite + オブジェクトストレージ分離

```
scores.sqlite       ← metadata (protein_hash, compound_hash, score, source, ...)
poses/{hash}.gz     ← pose blob を個別ファイル
```

- スコア検索が SQL で柔軟
- pose_blob は object storage (S3 / local FS) に分離 → HDF5 を使わない

**評価**:
| 観点 | 評価 |
|---|---|
| サイズ | ○ SQLite ~10-20GB + pose 136GB = 156 GB |
| 実装容易さ | × 大改修（依存追加、repository 書換） |
| アクセス速度 | ○ SQL インデックス活用可能 |
| Resume 差分検出 | ◎ SQL LEFT JOIN で瞬時 |
| 並列書込み | △ SQLite は単一 writer 制約 |

### 推奨: **案 A（Protein 束ね構造）**

**理由**:
1. サイズ目標達成（136 GB ≤ 200 GB）
2. 既存の `hdf5_docking_result_repository.py` への改修で済む（案 C のような大変更不要）
3. Phase 2 のメイン専有 write パターンを温存可能
4. Phase 4 で並列書込が必要になったら案 B へ拡張できる（段階的）

案 B は Phase 4 実運用で単一 HDF5 のファイルサイズがツール扱いにくい規模（100GB+）になった場合に追加適用する予備案。

---

## 3. 案 A の実装コスト見積もり

### 変更対象ファイル

| ファイル | 変更規模 | 内容 |
|---|---|---|
| `infrastructure/repositories/hdf5_docking_result_repository.py` | **大改修**（~150→200行） | save/load/find/has の全メソッドを protein-bundle 対応 |
| `tests/infrastructure/repositories/test_hdf5_docking_result_repository.py` | **大改修** | 全テストをスキーマ変更に追従 |
| `infrastructure/repositories/migrate_v2_to_v3.py` | **新規** | 既存 HDF5 → 新スキーマ変換スクリプト（optional） |
| `docking/screening_runner.py` | **軽微** | `_save_to_hdf5` / `_filter_unprocessed` を新 API に合わせる |

### 工数見積もり

| subtask | 内容 | 足軽 | 工数 |
|---|---|---|---|
| hdf5_repo_rewrite | repository 書換（save/load/has） | 1 | 1.5日 |
| hdf5_repo_tests | test 書換（~20件） | 1 | 1日 |
| screening_runner_adapt | Runner 側 API 追従 | 1 | 0.3日 |
| migration_script | 既存HDF5 migration（optional、殿判断） | 1 | 0.5日 |
| e2e_verification | 10×100 再実行で回帰確認 | 1 | 0.5日 |
| **合計** | — | **2足軽並列** | **実働 ~2日** |

### リスク

- HDF5 vlen bytes は h5py の bug で稀に問題になる（pose_blobs 配列に可変長バイト列を格納）。代替として固定サイズ (例: 1KB) の padding で対応可
- compound_hash set をメモリロード 6GB は devcontainer の RAM 制約（通常 8-16GB）で収まるが、本番インスタンスは 32GB+ 推奨

---

## 4. Phase 4 との整合性

- 殿ご裁可 Q1-A（A100/H100 クラウド）で 32-128GB RAM インスタンス利用可 → 差分 set 6GB は問題なし
- Dask worker 数 4-8 でも、メイン専有 write パターンなら案 A で十分。案 B は 16+ worker 規模で必要
- Phase 2 のスキーマ設計（source, top_n）は案 A でも維持（全 protein ダンプで計算可能）

---

## 5. 推奨アクション

### cmd_012 内で着手
- `subtask_012_hdf5_rewrite` を起票（上記 5 sub-subtask、2足軽並列 2日）
- Phase 4 本番着手の必須前提条件として位置付け

### 将来 Phase 4 運用時に再検討
- 単一 HDF5 ファイルが 100GB を超えたら案 B（UniProt shard）へ拡張
- ネット越し読み込みが多いなら案 C（SQLite + object storage）を検討

---

## 6. 結論

Phase 4 の HDF5 overhead 1339GB 問題は **案 A（Protein 束ね構造）** で **136GB（目標200GB 以下）** に削減可能。実装コストは 2 足軽並列で実働 2 日、既存 repository の書換が中心で設計変更は最小。cmd_012 Wave 2 で即着手推奨。
