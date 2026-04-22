# Phase 2 実装報告書

**作成日**: 2026-04-22
**対象**: docking_automation v2 Phase 2 (cmd_009)
**担当**: Track 009-F (足軽2号)

## 実装サマリ

| コンポーネント | 実装内容 | commit |
|---|---|---|
| ScreeningRunner 骨格 | N×M ドッキング司令塔、再開可能・冪等、Wave 1 シリアルスタブ | 4250f02 |
| HDF5 スキーマ更新 | pose_blob gzip-9圧縮、float32スコア、SWMR対応 | 2772320 |
| GridBoxCache missing_policy | skip/error/fallback_centroid の3モード対応 | 68bb6c2 |
| Dask統合 + JSONL | LocalCluster並列実行 + JSONL進捗ログ | 57ab931 |
| バグ修正 | MoleculeConverter に sdf_to_pdbqt を追加 (Track 009-F) | - |

## E2E 実行結果

| 項目 | 値 |
|---|---|
| タンパク質数 | 10 |
| 化合物数 | 10 |
| 総ペア数 | 100 |
| Run1 new_pairs | 100 |
| Run1 failed_pairs | 0 |
| Run2 new_pairs (resume確認) | 0 |
| Run2 reused_pairs (resume確認) | 100 |
| Run1 elapsed | 162.0 s |
| HDF5 サイズ | 702.1 KB (718,947 bytes) |
| Phase 4 外挿 (2e8ペア) | ~1,339 GB |

## タンパク質データ

- ソース: AlphaFold DB マウスプロテオーム (v6)
- ディレクトリ: `examples/input/afdb_mouse/`
- 選定: `protein_list.json` 先頭10件

| タンパク質ID | GridBox center (x, y, z) |
|---|---|
| AF-P10601-F1-model_v6 | (-4.35, -18.56, -2.40) |
| AF-P10889-F1-model_v6 | (-6.23, 1.69, 14.37) |
| AF-Q6W5C0-F1-model_v6 | (5.27, 7.49, -4.81) |
| AF-Q8JZM2-F1-model_v6 | (37.96, -7.91, -39.45) |
| AF-O35671-F1-model_v6 | (13.09, 10.07, -0.90) |
| AF-O54879-F1-model_v6 | (3.81, -9.19, -10.06) |
| AF-P0C0A3-F1-model_v6 | (-17.97, -14.46, 18.49) |
| AF-A2AF53-F1-model_v6 | (7.79, 3.45, -16.12) |
| AF-O08717-F1-model_v6 | (14.34, 8.68, -18.92) |
| AF-O35719-F1-model_v6 | (19.67, 7.40, -6.36) |

GridBox は fpocket (pocket_rank=1) で自動計算。一部タンパク質でサイズ最小値10Åに自動調整。

## 化合物データ

- ソース: `examples/input/ALDR/actives_subset.sdf` (10化合物)
- ALDR (aldose reductase) アクティブ化合物サブセット

## バグ修正: sdf_to_pdbqt 未実装

Track 009-C で実装された `dock_one_protein` が `MoleculeConverter.sdf_to_pdbqt()` を呼び出していたが、
該当メソッドが未実装だったため全100ペアが `compound_pdbqt_conversion_failed` エラーで失敗していた。

**修正内容** (`docking_automation/converters/molecule_converter.py`):
```python
def sdf_to_pdbqt(self, sdf_path: Path, output_path: Path) -> Path:
    """単一分子SDFファイルをPDBQTに変換する。"""
    suppl = Chem.SDMolSupplier(str(sdf_path), removeHs=False)
    mol = next((m for m in suppl if m is not None), None)
    if mol is None:
        raise ValueError(f"SDFファイルから有効な分子を読み込めません: {sdf_path}")
    return self._convert_to_pdbqt_with_meeko(mol, output_path)
```

## 計測値と考察

### スループット
- 100ペア / 162秒 ≒ 0.62ペア/秒 (Dask 4ワーカー、exhaustiveness=1)
- タンパク質ごとのDaskタスク分割により並列化

### HDF5ストレージ
- 7,189 bytes/ペア (gzip-9圧縮 pose_blob + メタデータ)
- Phase 4 スケール (2×10⁸ペア) では ~1,339 GB が必要
- 圧縮最適化や分散ストレージの検討が必要

### 冪等性・resume
- Run2 はHDF5の既存エントリを正確に検出し全100ペアをreuse
- 0.06秒で完了 (フィルタリングコストのみ)
- 冪等性が完全に機能していることを確認

## 残課題

- **ストレージ最適化**: Phase 4 (~1,339 GB) に向けた圧縮・シャーディング戦略が必要
- **スループット向上**: exhaustiveness=1 は検証用設定。本番では exhaustiveness=8〜32 が必要
- **グリッドボックス精度**: fpocket による自動計算。既知リガンドがある場合は結晶リガンド由来のグリッドを優先
- **エラーハンドリング強化**: JSONL ログに error フィールドを記録する改善が必要 (現状ログに error 詳細なし)
- **Phase 3**: 大規模スクリーニング (全マウスプロテオーム × 化合物ライブラリ) に向けた設計

## 生成ファイル

| ファイル | 説明 |
|---|---|
| `examples/phase2_e2e.py` | E2Eスクリプト |
| `examples/output/phase2_e2e.h5` | E2E実行結果 HDF5 (702 KB) |
| `examples/output/phase2_e2e_grid_cache.json` | fpocket GridBoxCache (10件) |
| `examples/output/phase2_e2e_screening.jsonl` | JSONL進捗ログ (200行) |
