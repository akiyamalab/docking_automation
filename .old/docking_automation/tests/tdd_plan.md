# docking_automation パッケージのためのTDD計画

このドキュメントは、docking_automationパッケージに対するテスト駆動開発(TDD)の計画を記述します。既存のコードベースと現状のテストケースを分析し、今後のテスト開発の方針を示します。

## 1. 現在のコード構造

### 1.1 パッケージ構成

docking_automationは以下の主要なモジュールで構成されています：

```
docking_automation/
├── application/         # アプリケーションサービスとワークフロー
│   ├── services/       # 結果解析などのサービス
│   └── workflows/      # ドッキング計算ワークフロー
├── docking/            # ドッキング計算関連のクラス
│   ├── entity/         # ドッキングエンティティ（タスク、結果など）
│   ├── service/        # ドッキングサービス
│   └── value_object/   # 値オブジェクト（設定、パラメータなど）
├── molecule/           # 分子構造関連のクラス
│   ├── entity/         # 分子エンティティ（化合物、タンパク質など）
│   ├── service/        # 分子準備サービス
│   └── value_object/   # 値オブジェクト（構造、形式など）
├── infrastructure/     # インフラストラクチャ
│   ├── compute/        # 計算リソース
│   ├── formats/        # ファイル形式変換
│   ├── persistence/    # データ永続化
│   └── tools/          # 外部ツール連携
└── presentation/       # プレゼンテーション層
    ├── api/            # API
    ├── cli/            # コマンドラインインターフェース
    └── web/            # Webインターフェース
```

### 1.2 主要クラス階層

#### 分子モジュール

- **Molecule** (抽象基底クラス)
  - 基本プロパティ: id, structure, format, properties, path, metadata
  - 抽象メソッド: prepare()
  - その他のメソッド: convert_to(), validate(), get_property()など

- **Compound** (化合物)
  - 基本プロパティ: id, path, structure, format, metadata
  - メソッド: add_metadata()
  
- **Protein** (タンパク質)
  - 基本プロパティ: id, path, structure, format, chains, metadata
  - メソッド: add_metadata(), has_active_site_defined()

#### ドッキングモジュール

- **DockingTask** (ドッキング計算タスク)
  - 基本プロパティ: id, ligand, receptor, configuration, status
  - メソッド: update_status(), is_ready(), prepare_for_execution()など
  
- **Ligand** (リガンド)
  - 基本プロパティ: compound
  
- **Receptor** (レセプター)
  - 基本プロパティ: protein

- **DockingResult** (ドッキング計算結果)
  - 基本プロパティ: id, task, poses, scores, created_at, execution_time

#### アプリケーションモジュール

- **SimpleDockingWorkflow** (シンプルなドッキング計算ワークフロー)
  - メソッド: execute() - 単一の化合物とタンパク質のドッキング計算
  
- **BatchDockingWorkflow** (バッチドッキング計算ワークフロー)
  - メソッド: execute() - 複数の化合物とタンパク質のドッキング計算

- **ResultAnalysisService** (結果解析サービス)
  - メソッド: calculate_score_statistics(), generate_score_histogram(), rank_results_by_score(), compare_results()

## 2. 既存のテストケース

### 2.1 テストディレクトリ構造

```
tests/
├── __init__.py
├── conftest.py        # テスト環境のセットアップとフィクスチャ
├── run_tests.sh       # テスト実行スクリプト
├── integration/       # 統合テスト
│   ├── __init__.py
│   ├── test_batch_workflow.py
│   ├── test_batch_workflow_mock.py
│   └── test_batch_workflow_run.py
├── scenario/          # シナリオテスト (現在空)
│   └── __init__.py
└── unit/              # 単体テスト
    ├── __init__.py
    ├── test_import.py
    ├── test_meeko.py
    ├── test_new_structure.py
    └── test_new_structure_direct.py
```

### 2.2 テストフィクスチャ (conftest.py)

現在のconftest.pyには以下のフィクスチャが定義されています：

- `setup_path()`: プロジェクトのルートディレクトリをパスに追加
- `temp_dir()`: テスト用の一時ディレクトリを提供
- `mock_compound()`: テスト用のモック化合物データを提供
- `mock_protein()`: テスト用のモックタンパク質データを提供

### 2.3 単体テスト

- **test_import.py**: パッケージのインポートと基本ディレクトリ構造の確認
- **test_meeko.py**: Meekoライブラリとの連携をテスト
- **test_new_structure.py**: 新しいディレクトリ構造への移行確認テスト
- **test_new_structure_direct.py**: 構造の直接的なテスト

### 2.4 統合テスト

- **test_batch_workflow.py**: バッチドッキングワークフローとその結果解析のテスト
- **test_batch_workflow_mock.py**: モックを使用したワークフローテスト
- **test_batch_workflow_run.py**: 実際にワークフローを実行するテスト

### 2.5 シナリオテスト

現在は空のディレクトリで、特定のユースケースを検証するテストはまだ実装されていないようです。

## 3. テスト状況の評価

### 3.1 現状の課題

1. **多くのクラスに対するユニットテストの不足**
   - エンティティクラスの個別テストが少ない
   - 値オブジェクトの個別テストが少ない

2. **テスト範囲の制限**
   - 一部の機能のみがテストされている
   - エラーケースやエッジケースのテストが不足

3. **テストカバレッジの偏り**
   - 既存のテストはワークフローの一部のみをカバー
   - インフラストラクチャ層のテストが少ない

## 4. TDDのための提案

### 4.1 追加すべきテストケース

#### 4.1.1 エンティティクラスのテスト

**Moleculeクラス**
- 初期化と基本プロパティのテスト
- validate()メソッドの検証
- get_property()とadd_property()のテスト
- エラーケースのテスト（無効な構造や形式）

**Compoundクラス**
- 初期化と名前の自動設定テスト
- add_metadata()の動作検証
- パスなしでの初期化テスト

**Proteinクラス**
- 初期化とチェーン情報のテスト
- has_active_site_defined()の検証
- 複数チェーンの処理テスト

**DockingTaskクラス**
- 初期化とIDの自動生成テスト
- ステータス変更テスト（updateメソッド）
- is_ready()の条件テスト
- タイムアウト処理のテスト
- マルチレセプター対応のテスト

**DockingResultクラス**
- 初期化と基本プロパティのテスト
- get_best_score()と他のスコア関連メソッドのテスト
- ポーズの取得と検証テスト

#### 4.1.2 値オブジェクトのテスト

**GridBoxクラス**
- 初期化と座標値のバリデーション
- サイズの計算と検証
- 無効な値のエラーハンドリング

**DockingConfigurationクラス**
- 初期化と基本プロパティのテスト
- validate()メソッドのテスト
- with_grid_box()メソッドの検証

**DockingParametersクラス**
- パラメータの追加と取得テスト
- to_dict(), from_dict()の検証
- パラメータの検証ロジックテスト

**Scoreクラス**
- 初期化と基本プロパティのテスト
- 異なるスコアタイプの比較テスト

#### 4.1.3 サービスクラスのテスト

**MoleculePreparationServiceクラス**
- prepare_ligand()とprepare_receptor()のテスト
- 異なるツール用の準備テスト
- エラーケースのテスト

**MoleculePreparationFactoryクラス**
- create_for_tool()メソッドのテスト
- サポートされていないツールのエラーハンドリング

**DockingServiceクラス**
- create_task()のテスト
- execute()の処理テスト
- get_default_configuration()の検証

**ResultAnalysisServiceクラス**
- calculate_score_statistics()のテスト
- generate_score_histogram()の検証
- rank_results_by_score()と結果ランキングのテスト
- compare_results()と複数結果の比較テスト

#### 4.1.4 統合テスト

**SimpleDockingWorkflowクラス**
- execute()メソッドの全体フロー検証
- 異なる入力パターンのテスト
- エラーケースの処理検証

**BatchDockingWorkflowクラス**
- 複数の化合物と単一タンパク質のバッチ処理テスト
- 並列処理の動作確認
- 結果集約機能のテスト
- エラー処理とリカバリーテスト

#### 4.1.5 シナリオテスト

**基本的なドッキングシナリオ**
- 単一の化合物とタンパク質のドッキング計算
- 結果の検証と可視化

**複数リガンドドッキングシナリオ**
- 複数の化合物ライブラリを使用したドッキング
- ランキングと結果比較

**インフラストラクチャ連携シナリオ**
- 外部ツール（Vina, Meekoなど）との連携テスト
- 計算環境（ローカル、クラウドなど）の切り替えテスト

**CLIからの実行シナリオ**
- コマンドラインからのワークフロー実行テスト
- 引数処理と結果出力のテスト

### 4.2 TDD実施計画

#### フェーズ1: 基本的な単体テストの実装

1. 値オブジェクトのテスト実装
   - GridBox, DockingConfiguration, DockingParameters, Score

2. エンティティのテスト実装
   - Compound, Protein, Ligand, Receptor
   - DockingTask, DockingResult

3. 基本サービスのテスト実装
   - MoleculePreparationService
   - DockingService
   - ResultAnalysisService

#### フェーズ2: 高度な単体テストと統合テスト

1. エラーケースとエッジケースのテスト
   - 異常入力値の処理
   - タイムアウト処理
   - 並行処理のエッジケース

2. ワークフローの統合テスト
   - SimpleDockingWorkflow
   - BatchDockingWorkflow

3. サービス間連携のテスト
   - 準備サービスとドッキングサービスの連携
   - 結果解析サービスとの連携

#### フェーズ3: シナリオテストと検証

1. E2Eシナリオの実装
   - 実際のユースケースを再現するテスト

2. パフォーマンステスト
   - 大規模データセットでの処理検証

3. インフラストラクチャ連携テスト
   - 外部ツールとの連携

### 4.3 テスト優先順位

実装の優先順位としては、以下のように提案します：

1. 基本的なエンティティと値オブジェクトのテスト
   - Compound, Protein, GridBox, DockingConfiguration

2. 核となるサービスのテスト
   - MoleculePreparationService, DockingService

3. ワークフローの統合テスト
   - SimpleDockingWorkflow, BatchDockingWorkflow

4. エラーケースとエッジケースのテスト
   - タイムアウト、無効な入力、境界値など

## 5. テスト実装例

### 5.1 Compoundクラスのテスト例

```python
import pytest
from docking_automation.molecule.entity.compound import Compound

def test_compound_initialization():
    """Compoundの初期化をテスト"""
    # 基本初期化
    compound = Compound(
        id="test_id",
        path="path/to/compound.sdf"
    )
    
    assert compound.id == "test_id"
    assert compound.path == "path/to/compound.sdf"
    assert compound.name == "compound"  # ファイル名からの自動設定
    assert compound.metadata == {}
    
def test_compound_custom_name():
    """名前のカスタム設定をテスト"""
    compound = Compound(
        id="test_id",
        path="path/to/compound.sdf"
    )
    
    # 名前の変更
    compound.name = "custom_name"
    assert compound.name == "custom_name"
    
def test_compound_add_metadata():
    """メタデータ追加をテスト"""
    compound = Compound(
        id="test_id",
        path="path/to/compound.sdf"
    )
    
    # メタデータ追加
    compound.add_metadata({"molecular_weight": 180.2})
    assert "molecular_weight" in compound.metadata
    assert compound.metadata["molecular_weight"] == 180.2
    
    # 追加のメタデータ
    compound.add_metadata({"logP": 2.5})
    assert "logP" in compound.metadata
    assert "molecular_weight" in compound.metadata  # 既存のメタデータも保持
```

### 5.2 GridBoxクラスのテスト例

```python
import pytest
from docking_automation.docking.value_object.grid_box import GridBox

def test_grid_box_initialization():
    """GridBoxの初期化をテスト"""
    grid_box = GridBox(
        center_x=10.0,
        center_y=20.0,
        center_z=30.0,
        size_x=15.0,
        size_y=25.0,
        size_z=35.0
    )
    
    assert grid_box.center_x == 10.0
    assert grid_box.center_y == 20.0
    assert grid_box.center_z == 30.0
    assert grid_box.size_x == 15.0
    assert grid_box.size_y == 25.0
    assert grid_box.size_z == 35.0
    
def test_grid_box_bounds():
    """GridBoxの境界計算をテスト"""
    grid_box = GridBox(
        center_x=10.0,
        center_y=20.0,
        center_z=30.0,
        size_x=10.0,
        size_y=20.0,
        size_z=30.0
    )
    
    min_x, max_x = grid_box.get_x_bounds()
    min_y, max_y = grid_box.get_y_bounds()
    min_z, max_z = grid_box.get_z_bounds()
    
    assert min_x == 5.0  # center - size/2
    assert max_x == 15.0  # center + size/2
    assert min_y == 10.0
    assert max_y == 30.0
    assert min_z == 15.0
    assert max_z == 45.0
    
def test_grid_box_validation():
    """GridBoxのバリデーションをテスト"""
    # 有効なGridBox
    valid_grid = GridBox(
        center_x=10.0,
        center_y=20.0,
        center_z=30.0,
        size_x=15.0,
        size_y=25.0,
        size_z=35.0
    )
    assert valid_grid.validate() == True
    
    # 無効なGridBox（サイズが負）
    invalid_grid = GridBox(
        center_x=10.0,
        center_y=20.0,
        center_z=30.0,
        size_x=-5.0,  # 負のサイズ
        size_y=25.0,
        size_z=35.0
    )
    assert invalid_grid.validate() == False
```

## 6. まとめ

このTDD計画では、docking_automationパッケージの現状を分析し、テストカバレッジを向上させるための具体的な方針を示しました。段階的にテストを実装し、コードの品質と信頼性を向上させていくことが重要です。

優先順位としては、基本的なエンティティと値オブジェクトのテストから始め、徐々に複雑なサービスとワークフローのテストへと進むことを推奨します。また、エラーケースとエッジケースのテストも忘れずに実装することが重要です。

TDDの基本的なサイクルを守り、テストの作成、実装、リファクタリングを繰り返すことで、高品質なコードを維持しながら機能を拡張していくことができます。