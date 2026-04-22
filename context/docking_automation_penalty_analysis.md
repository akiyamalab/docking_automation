# Uni-Dock Penalty Score 根本原因分析

**作成日**: 2026-04-22
**対象**: docking_automation Phase 3 → Phase 4 準備（cmd_012 subtask_012_gunshi_a タスクA）
**作成**: 軍師

---

## 1. 観察データ

### 1.1 Phase 3 再実行結果（commit 7d99d41）

| 区分 | 件数 | 割合 |
|---|---|---|
| 新規成功（valid score） | 690 | 69.0% |
| 失敗（penalty filter + output_missing） | 310 | 31.0% |
| 計 | 1000 | 100% |

内訳（足軽3号報告より推定）:
- `unidock_output_missing`（out_pdbqt ファイル未生成）
- `unidock_score_filtered`（score ≥ 5.0 または ≤ -30.0）

### 1.2 スコア分布（フィルタ前の観測値）

| 統計 | 値 (kcal/mol) |
|---|---|
| 有効 690 件 mean | -3.90 |
| 有効 690 件 min | -7.31 |
| 有効 690 件 max | 4.87 |
| **異常値例**（フィルタ対象） | 1,623,719（指数的ペナルティ）、大量の正値 |

Vina でも同データで 31/100 件が正値（≥0）→ **ペアそのものがドッキング困難**（リガンド-グリッド不整合）が主因の可能性が高い。

---

## 2. 失敗パターンの仮説

### 仮説 H1: グリッドボックスが小さすぎる

- `GridBox.from_fpocket` + `fpocket` が返すポケット領域は典型的に 15-25 Å 立方
- リガンド最大径がグリッド内に収まらない場合、Uni-Dock は pose 生成不能 → penalty
- ALDR actives_final のリガンドは MW ~300-600 で分子長 10-20 Å あり、境界ケース

**検証方法**:
- 失敗ペアのリガンド分子長（radius of gyration × 2）と GridBox size を比較
- `size < 2 × ligand_radius + 4Å` のペアで失敗率が高いか

**対策候補**:
- GridBox size に padding = 5-8 Å を追加して 20-30 Å 立方を標準化
- ligand サイズ依存の動的 padding

### 仮説 H2: PDBQT 変換の部分破損

- Meeko `MoleculePreparation` は特殊官能基（マクロサイクル、リン酸、硫酸エステル）で警告を出すことがある
- 変換失敗を気付かず PDBQT が空/不完全になり、UniDock 側で score が出力されない
- 足軽3号の報告には「unidock_output_missing」が含まれており、これが該当する疑い

**検証方法**:
- 失敗 PDBQT リガンドの原子数・回転可能結合数を確認
- Meeko warning をキャプチャ（stderr）、失敗ペアと相関

**対策候補**:
- compound_pipeline の Meeko 前処理で warning を error 扱いに厳格化
- PDBQT サイズが閾値未満なら前処理段階で failed 記録

### 仮説 H3: UniDock search_mode=balance での精度限界

- `search_mode='balance'` は速度重視（exhaustiveness 相当の途中値）
- 柔軟リガンド（RotB ≥ 10）は balance では収束不足で penalty 返却の可能性
- 論文（JCTC 2023）では detail モードで同じリガンドを再実行すると rescue される例あり

**検証方法**:
- 失敗ペアのうちランダム 20 件を `search_mode='detail'` で再実行
- rescue 率を測定

**対策候補**:
- 2段階実行: balance で 1回目 → penalty 返却分だけ detail で 2回目
- コストは 1.5-2× だが失敗率を半減できる見込み

### 仮説 H4: リガンド自体がドッキング不能（生物学的）

- ALDR actives_final は ALDR 阻害剤として設計されたが、他タンパク質との組合せでは相互作用しない可能性
- N×M スクリーニングの本質: 大半は「結合しない」ペア → 物理的に pose が得られない
- この場合 penalty は "正しい失敗" であり、修正対象ではない

**検証方法**:
- 失敗ペアの Vina 側のスコアも同様に positive なら H4
- Vina が -6〜-8 を返すのに UniDock が penalty を返すなら UniDock 固有の問題

**対策候補**:
- 原則修正なし。ただし `failed_pairs` を `no_binding` と `docking_error` に分類する分析レイヤを追加
- Phase 4 では penalty = no_binding として扱い、scoring distribution 分析に含める

---

## 3. 仮説の優先度と推定寄与率

| 仮説 | 推定寄与率 | 根拠 |
|---|---|---|
| H4 (no_binding) | 40-60% | Vina も類似に positive score を返している事実。マウス全 protein × 薬物様 compound は大半が非結合 |
| H1 (grid size) | 15-25% | ALDR リガンドは中サイズ、fpocket デフォルト grid で境界 |
| H2 (PDBQT破損) | 5-15% | Meeko warning の頻度に依存、要実測 |
| H3 (search_mode balance 限界) | 10-20% | 柔軟リガンドで顕在化、10-20% 程度と推定 |

---

## 4. 推奨対策（Phase 4 着手前）

### 短期（cmd_012 内で実施可）

1. **対策-1: グリッド padding +5Å 実験**
   - GridBoxCache の GridBox 生成時に size をデフォルト +5Å 拡張
   - 失敗ペア 310 件のうち H1 起因分を rescue
   - コスト: GridBoxCache JSON 再生成のみ、実装工数 小

2. **対策-2: search_mode 2段階実行**
   - ScreeningRunner に `rescue_mode` オプション追加: penalty ペアのみ `detail` で再実行
   - Phase 4 本番の `failed_pairs` の subset に適用
   - コスト: 実装 0.5日、計算時間 +10-20%

3. **対策-3: Meeko warning 厳格化**
   - compound_pipeline で warning → error 昇格、前処理段階で除外
   - H2 起因分を事前に除去（HDF5 cache に `no_binding` と `prep_failed` を明示区分）

### 中長期（Phase 4 本番運用）

4. **failed_pairs 分類強化**
   - `failed_reason` を enum 化: `no_binding`, `grid_too_small`, `prep_failed`, `timeout`, `unknown`
   - HDF5 metadata カラムに記録、後段の統計解析で区別

5. **Vina との 2-backend 検証**
   - cmd_011 で証明済の r=0.9580 を活用: UniDock 失敗ペアを Vina で rescue
   - 計算コスト vs 回収率のトレードオフ要検討

---

## 5. Phase 4 での影響見積もり

### 現状ベースライン（対策なし）
- 失敗率 31% → 2×10⁸ ペアで **6.2×10⁷ 件** が penalty
- 有効データ 1.38×10⁸ 件 → 十分な analysis セット

### 対策-1/2 適用時（推定）
- 失敗率 15-20% に改善見込み
- 有効データ 1.6-1.7×10⁸ 件 → analysis 品質向上

### リスク
- 対策を打たずに Phase 4 を走らせると、後から「失敗ペアを救済したい」場合に 6.2×10⁷ 件の再計算となり、GPU コスト大幅増
- 一方、対策-1 (padding) は下手に大きくするとドッキング時間が増える（grid サイズが体積比で効く）ため最小限に留めるべき

---

## 6. 結論

- penalty 発生の 40-60% は生物学的 non-binding（修正不要）と推定
- 残り 40-60% のうち grid size / search_mode / Meeko 前処理の改善で半減可能
- **Phase 4 着手前に対策-1（grid padding +5Å）と対策-2（2段階実行）の実装を推奨**
- `failed_reason` 分類を HDF5 に追加しておけば Phase 4 後でも補助分析で原因特定できる

## 7. 次タスク候補

- `subtask_012_x_penalty_investigate`: 失敗 310 件のリガンド分子量・RotB・grid 関係実測（足軽1名 0.5日）
- `subtask_012_y_grid_padding`: GridBoxCache.from_fpocket に padding オプション追加（足軽1名 0.3日）
- `subtask_012_z_rescue_mode`: ScreeningRunner rescue_mode 実装（足軽1名 0.7日）
