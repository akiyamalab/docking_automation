# Docking Automation

タンパク質-化合物ドッキング計算の自動化システム（ドメイン駆動設計アプローチ）

## 概要

このプロジェクトは、タンパク質-化合物ドッキング計算を効率的に自動化するためのシステムです。ドメイン駆動設計（DDD）の原則に基づいて設計されており、以下の特徴を持ちます：

- 複数のドッキングツール（AutoDock Vina, gnina, diffdock, REstretto等）のサポート
- 柔軟なデータ保存形式（CSV, JSON, RDB, NoSQL）
- 大規模計算のための分散処理対応（ローカル実行、SLURM、クラウド）
- 明確なドメインモデルによる拡張性と保守性の高さ

## 主要機能

- 分子準備：化合物とタンパク質の前処理
- ドッキング実行：様々なドッキングツールによる計算
- ジョブ管理：計算ジョブの分散と監視
- 結果管理：計算結果の解析と評価
- データ管理：複数形式での永続化と取得

## インストール

```bash
# 開発版としてインストール
git clone https://github.com/yourusername/docking_automation.git
cd docking_automation
pip install -e .

# または直接インストール
pip install docking_automation
```

## 必要条件

- Python 3.8以上
- 以下の外部ツール（別途インストールが必要）：
  - AutoDock Vina
  - Open Babel
  - MGLTools (prepare_ligand, prepare_receptor)

## 使用例

### 単一ドッキング計算

```python
from docking_automation.application.workflows.simple_docking import SimpleDockingWorkflow

# ワークフローの実行
workflow = SimpleDockingWorkflow()
result = workflow.execute(
    ligand_path="path/to/ligand.sdf",
    receptor_path="path/to/receptor.pdb",
    config={
        "center_x": -26.75,
        "center_y": 5.82,
        "center_z": -8.13,
        "size_x": 20.0,
        "size_y": 20.0,
        "size_z": 20.0,
        "exhaustiveness": 8,
    }
)

# 結果の処理
print(f"Best score: {result.scores[0].value}")
result.save("output/result.csv")
```

### バッチ処理

```python
from docking_automation.application.workflows.batch_docking import BatchDockingWorkflow

# バッチワークフローの実行
workflow = BatchDockingWorkflow()
results = workflow.execute(
    ligands_dir="path/to/ligands",
    receptor_path="path/to/receptor.pdb",
    output_dir="output",
    n_jobs=4  # 並列数
)

# 結果のランキング
ranked_results = results.rank_by_score()
ranked_results.save("output/ranked_results.csv")
```

## ライセンス

MIT

## 貢献

プルリクエストは歓迎します。大きな変更を加える場合は、まずissueを開いて変更内容を議論してください。

## 開発者

- あなたの名前 – [@yourusername](https://github.com/yourusername)