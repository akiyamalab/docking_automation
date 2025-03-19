# Docking Automation Framework

分子ドッキング計算を自動化するためのフレームワーク

## 概要

Docking Automation Frameworkは、タンパク質-化合物ドッキング計算を簡単かつ効率的に実行するためのツールです。創薬研究や構造生物学研究において、複数の化合物をタンパク質に対してドッキングシミュレーションする作業を自動化し、再現性と効率性を高めます。

### 主な特徴

- **簡単な操作**: 共通インターフェースによる複雑なドッキング計算の実行
- **高い再現性**: 明確に定義されたワークフローによる計算の完全な再現性
- **効率的な処理**: 並列処理による大規模計算の高速化
- **複数ツール対応**: 様々なドッキングツールを統一インターフェースで利用可能

## システム要件

- **Python**: 3.8以上
- **対応OS**: Linux, macOS, Windows
- **必要なディスク容量**: 500MB以上（ドッキング結果保存用の追加容量が必要）
- **推奨メモリ**: 8GB以上（大規模計算には16GB以上推奨）

## インストール方法

### 基本インストール

```bash
# 基本パッケージのインストール
pip install -e .
```

### 追加機能のインストール

```bash
# 化学情報学ライブラリを含めたインストール
pip install -e .[chem]

# ドッキング関連ライブラリを含めたインストール
pip install -e .[docking]

# すべての機能を含めたインストール
pip install -e .[full]
```

### RDKitのインストール

RDKitはconda経由でのインストールが推奨されています：

```bash
conda install -c conda-forge rdkit
```

## 基本的な使用方法

### 主要なクラスと機能

- **Protein**: タンパク質構造を表現するクラス
- **CompoundSet**: 複数の化合物をまとめて扱うクラス
- **GridBox**: ドッキング計算の探索空間を定義するクラス
- **AutoDockVina**: AutoDock Vinaを使用したドッキング計算を実行するクラス
- **DockingResult**: ドッキング計算の結果を管理するクラス

### 基本的なワークフロー

1. **タンパク質と化合物の準備**: 
   - タンパク質構造ファイル（PDB, MOL2など）の読み込み
   - 化合物ファイル（SDF, MOL2など）の読み込み

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

以下は、基本的なドッキング計算を実行する例です：

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

# 結果の表示
for i, result in enumerate(top_hits):
    print(f"{i+1}. Score: {result.docking_score}, Compound: {result.compound_id}")
    print(f"   Pose file: {result.result_path}")
```

より詳細な例は[examples/simple_docking.py](examples/simple_docking.py)を参照してください。

## 応用例

### 異なるドッキングツールの使用

```python
# AutoDock Vinaの代わりにREStrettoを使用
from docking_automation.docking import REStrettoDocking

docking_tool = REStrettoDocking()
results = docking_tool.run_docking(protein, compounds, grid_box)
```

### バッチ処理と並列計算

```python
from docking_automation.infrastructure.executor import ParallelExecutor

# 並列実行エンジンの設定
executor = ParallelExecutor(max_workers=4)

# バッチ処理の実行
results = docking_tool.run_docking(
    protein, 
    compounds, 
    grid_box, 
    executor=executor,
    batch_size=100  # 100化合物ごとにバッチ処理
)
```

### 結果の保存と読み込み

```python
# 結果をCSVファイルに保存
results.save_to_csv("docking_results.csv")

# 結果をJSONファイルに保存
results.save_to_json("docking_results.json")

# 結果の読み込み
from docking_automation.docking import DockingResultCollection
loaded_results = DockingResultCollection.load_from_csv("docking_results.csv")
```

## トラブルシューティング

### インストール時の問題

#### RDKitのインストールエラー

RDKitのインストールに問題がある場合は、Condaを使用したインストールをお試しください：

```bash
conda install -c conda-forge rdkit
pip install -e .
```

#### OpenBabelのインストールエラー

OpenBabelのインストールに問題がある場合：

```bash
# Linuxの場合
sudo apt-get install openbabel
# または
conda install -c conda-forge openbabel

# macOSの場合
brew install open-babel
```

### 実行時のエラー

#### メモリエラー

大規模な化合物セットを処理する際にメモリエラーが発生する場合は、バッチサイズを小さくして実行してください：

```python
results = docking_tool.run_docking(protein, compounds, grid_box, batch_size=50)
```

#### タンパク質準備のエラー

タンパク質の前処理中にエラーが発生する場合は、以下を確認してください：

- PDBファイルの形式が正しいか
- 欠損残基や原子がないか
- 水素原子の追加が必要か

```python
# 水素原子を明示的に追加
protein = Protein(protein_path, add_hydrogens=True)
```

## サポートとリソース

- **ドキュメント**: [完全なドキュメント](https://docking-automation.readthedocs.io/)
- **バグ報告**: [Issue Tracker](https://github.com/yourusername/docking_automation/issues)
- **ディスカッション**: [GitHub Discussions](https://github.com/yourusername/docking_automation/discussions)
- **連絡先**: example@email.com

## ライセンス

このプロジェクトはMITライセンスの下で公開されています。詳細は[LICENSE](LICENSE)ファイルをご覧ください。