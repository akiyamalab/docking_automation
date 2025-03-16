# Docking Automation Framework

分子ドッキング計算を自動化するためのドメイン駆動設計フレームワーク

## 実装状況

現在、以下の実装が完了しています：

1. **基本クラス構造の実装**：
   - 新規作成: docking_automation/converters/molecule_converter.py
   - 新規作成: docking_automation/infrastructure/executor/task.py
   - 新規作成: docking_automation/infrastructure/executor/task_manager.py
   - 新規作成: docking_automation/docking/docking_result_collection.py
   - 修正: docking_automation/molecule/protein.py
   - 修正: docking_automation/molecule/compound_set.py
   - 修正: docking_automation/docking/docking_result.py
   - 修正: docking_automation/docking/grid_box.py
   - 修正: docking_automation/infrastructure/repositories/docking_result_repository.py
   - 修正: docking_automation/infrastructure/executor/executor.py
   - 修正: docking_automation/infrastructure/executor/sequential_executor.py
   - 修正: docking_automation/infrastructure/executor/slurm_executor.py

2. **テストフレームワークの実装**：
   - tests/molecule/test_protein.py
   - tests/molecule/test_compound_set.py
   - tests/converters/test_molecule_converter.py
   - tests/docking/test_grid_box.py
   - tests/docking/test_docking_result.py
   - tests/docking/test_docking_result_collection.py
   - tests/infrastructure/repositories/test_docking_result_repository.py
   - tests/infrastructure/executor/test_task.py
   - tests/infrastructure/executor/test_task_manager.py
   - tests/infrastructure/executor/test_sequential_executor.py
   - tests/infrastructure/executor/test_slurm_executor.py
   - tests/integration/test_docking_workflow.py
   - tests/conftest.py

現時点では、すべてのメソッドはNotImplementedErrorを発生させるようになっており、具体的な実装は行っていません。これにより、設計の骨組みとテストフレームワークが整い、今後の実装作業の基盤が整いました。

詳細な設計と実装状況については、[docking_automation_design_improvements.md](docking_automation_design_improvements.md)を参照してください。

## インストール方法

### 依存関係

本プロジェクトは以下の依存関係を持ちます：

- Python 3.8以上
- RDKit: 化学情報学ライブラリ
- OpenBabel: 分子ファイル形式変換
- meeko: リガンド準備（AutoDock Vina用）
- sortedcontainers: ソート済みコンテナ
- pytest: テストフレームワーク

### インストール手順

開発モードでインストールするには、以下のコマンドを実行します：

```bash
# 基本パッケージのインストール
pip install -e .

# 開発用依存関係を含めたインストール
pip install -e .[dev]

# 化学情報学ライブラリを含めたインストール
pip install -e .[chem]

# すべての依存関係を含めたインストール
pip install -e .[full]
```

## 開発環境のセットアップ

1. リポジトリをクローン：
   ```bash
   git clone https://github.com/yourusername/docking_automation.git
   cd docking_automation
   ```

2. 依存関係のインストール：
   ```bash
   pip install -e .[full,dev]
   ```

3. テストの実行：
   ```bash
   pytest
   ```

## テストの実行

テストを実行するには、以下のコマンドを使用します：

```bash
# すべてのテストを実行
pytest

# 特定のテストファイルを実行
pytest tests/molecule/test_protein.py

# カバレッジレポートを生成
pytest --cov=docking_automation
```

## 1. プロジェクトビジョン

### 1.1 目的と背景

Docking Automation Frameworkは、タンパク質-化合物ドッキング計算の自動化プロセスを標準化・効率化することを目的としています。創薬研究の現場では、複数のドッキングツールを用いた大規模スクリーニングが一般的ですが、これらの作業フローはしばしば研究者ごとに異なり、再現性や効率性に課題があります。本フレームワークは以下の課題解決を目指します：

- **再現性の確保**: 明確に定義されたワークフローによる計算の完全な再現性
- **効率化**: スケーラブルな分散処理による大規模計算の高速化
- **標準化**: 様々なドッキングツールの統一インターフェース提供
- **知識の共有**: ドメインモデルによる計算化学・構造生物学の知識の明示的表現

### 1.2 想定利用者

- **計算化学者**: 大規模スクリーニングや複雑なドッキング戦略の構築
- **構造生物学者**: タンパク質構造に基づくリガンド相互作用解析
- **創薬研究者**: 新規化合物のスクリーニングと最適化
- **バイオインフォマティシャン**: 計算パイプラインへの組み込み

### 1.3 長期ビジョン

将来的には、以下の機能を提供するプラットフォームへと発展させることを目指します：

- AI支援による自動ドッキングパラメータ最適化
- 分子力学計算と統合した高精度相互作用解析
- クラウドネイティブな大規模分散処理
- 創薬プロセス全体を支援する包括的なワークフロー環境

## 2. アーキテクチャ

### 2.1 ドメイン駆動設計の採用理由

本プロジェクトは、以下の理由からドメイン駆動設計（DDD）アプローチを採用しています：

- **ドメイン複雑性の管理**: 分子ドッキングの複雑な概念と関係性を適切にモデル化
- **分野知識の明示化**: 暗黙知となりがちな専門知識をコードで明示的に表現
- **変更への適応性**: 技術的実装と独立したドメインモデルによる長期的保守性
- **拡張性**: 新たなドッキングアルゴリズムや分子表現を容易に追加できる構造

### 2.2 システムアーキテクチャ

本システムは、以下の4つの主要レイヤーで構成されるヘキサゴナルアーキテクチャを採用しています：

#### 2.2.1 ドメイン層
- システムの中心となる概念モデルと業務ルール
- 外部システムや技術的詳細に依存しない純粋なドメインロジック
- 主要概念：Protein, Ligand, DockingResult, BindingSite など

#### 2.2.2 アプリケーション層
- ユースケースを実現するためのワークフローとサービス
- ドメイン層を組み合わせた高レベル操作
- 主要概念：DockingWorkflow, BatchProcessor, ResultAnalyzer など

#### 2.2.3 インフラストラクチャ層
- 技術的詳細の実装
- 外部サービス・ツールとの連携
- ストレージ、分散処理、外部APIなど
- 主要概念：DockingToolAdapter, ResultRepository, ClusterExecutor など

#### 2.2.4 インターフェース層
- ユーザーとシステムの相互作用
- CLI, API, GUIなど
- 主要概念：CLICommands, RESTEndpoints など

### 2.3 モジュール構成

```
docking_automation/
├── domain/
│   ├── molecule/            # 分子関連の基本モデル
│   ├── docking/             # ドッキング計算のコアロジック
│   └── analysis/            # 結果解析のドメインロジック
├── application/
│   ├── services/            # アプリケーションサービス
│   ├── workflows/           # ワークフロー定義
│   └── dto/                 # データ転送オブジェクト
├── infrastructure/
│   ├── repositories/        # データ永続化実装
│   ├── executor/            # 分散実行エンジン
│   ├── tools/               # 外部ツール連携
│   └── adapters/            # 外部システムアダプタ
└── interface/
    ├── cli/                 # コマンドラインインターフェース
    ├── api/                 # プログラミングインターフェース
    └── web/                 # Webインターフェース（将来）
```

## 3. ドメインモデル設計

### 3.1 主要ドメインオブジェクト

本システムでは、以下の主要なドメインオブジェクトを定義します：

#### 3.1.1 分子モデル
- **Molecule**: 全ての分子の基底抽象クラス
  - **Protein**: 受容体となるタンパク質を表現
  - **Ligand**: 低分子化合物（リガンド）を表現
  - **CompoundSet**: 複数のリガンドをまとめたコレクション

#### 3.1.2 ドッキングモデル
- **BindingSite**: タンパク質の結合部位を表現
- **GridBox**: ドッキング計算の探索空間
- **DockingParameters**: ドッキング計算のパラメータ
- **DockingResult**: ドッキング計算の結果
  - **Pose**: ドッキングポーズ（配座）
  - **InteractionProfile**: 分子間相互作用のプロファイル
  - **Score**: スコアとその詳細内訳

#### 3.1.3 ワークフローモデル
- **WorkflowStep**: ワークフローの各ステップを表現
- **WorkflowState**: ワークフローの状態を管理
- **ExecutionContext**: 実行コンテキスト（環境設定など）

### 3.2 値オブジェクトとエンティティの区別

DDDでは、オブジェクトをその同一性に基づいて「エンティティ」と「値オブジェクト」に分類します：

#### エンティティ（ID で識別）
- Protein: PDB ID や内部IDで識別
- Ligand: SMILES, InChI, または内部IDで識別
- DockingJob: 一意のジョブIDで識別

#### 値オブジェクト（属性の値で識別）
- AtomCoordinate: x, y, z 座標
- GridBox: 中心座標とサイズ
- DockingParameters: 各種パラメータの集合

### 3.3 集約とバウンダリ

複雑さを管理するため、関連するオブジェクトをグループ化した「集約」を定義します：

- **Protein集約**: Protein, Chain, Residue, Atom
- **Ligand集約**: Ligand, Atom, Bond, FunctionalGroup
- **DockingResult集約**: DockingResult, Pose, Score, InteractionProfile

各集約のルートエンティティ（Protein, Ligand, DockingResultなど）を通じてのみ、集約内の他のオブジェクトにアクセスします。

## 4. 実装戦略

### 4.1 分子表現の拡張方針

現在の最小実装から以下の方向に拡張します：

- **内部表現の強化**:
  - RDKitベースの分子内部表現
  - 原子・結合情報の完全サポート
  - 分子特性計算機能の統合

- **インターコンバーション**:
  - 様々な分子ファイル形式間の変換
  - 内部表現と外部ツール間のデータ変換

- **検証メカニズム**:
  - 分子構造の妥当性検証
  - タンパク質準備状態の検証

### 4.2 ドッキングワークフローの設計

ドッキングプロセスを以下のフェーズに分割し、各フェーズをカスタマイズ可能にします：

1. **準備フェーズ**:
   - 受容体準備（水素付加、電荷割り当て、構造最適化）
   - リガンド準備（3D構造生成、コンフォメーション生成）
   - グリッドボックス定義（自動／手動）

2. **ドッキングフェーズ**:
   - 単一ツールによるドッキング
   - 複数ツールによる並列ドッキング
   - アンサンブルドッキング戦略

3. **後処理フェーズ**:
   - ポーズのクラスタリング
   - 分子間相互作用解析
   - 結果の可視化準備
   - 分子動力学シミュレーション準備

4. **評価フェーズ**:
   - スコアリングとランキング
   - ファーマコフォアベース評価
   - カスタムフィルタリング

各フェーズは「戦略パターン」を採用し、異なる実装を容易に切り替えられるようにします。

### 4.3 拡張性とプラグイン機構

システムの拡張性を高めるために、以下のプラグイン機構を実装します：

- **ドッキングツールプラグイン**:
  - 共通インターフェースを実装する新規ドッキングツールアダプタ
  - ツール固有のパラメータをサポートする柔軟な設計

- **スコア関数プラグイン**:
  - カスタムスコアリング関数の登録メカニズム
  - 複数スコアの組み合わせ機能

- **分析プラグイン**:
  - 結果解析用のカスタムアルゴリズム追加
  - 可視化プラグイン

プラグインの登録には、依存性注入（DI）コンテナと工場パターンを活用します。

### 4.4 実行エンジンの設計

複数の実行環境をサポートする柔軟な実行エンジンを実装します：

- **逐次実行（SequentialExecutor）**:
  - 単一マシンでの逐次処理
  - デバッグと小規模計算向け

- **ローカル並列（LocalParallelExecutor）**:
  - マルチプロセス/マルチスレッドによる並列処理
  - 単一マシンでの効率的な実行

- **クラスター実行（ClusterExecutor）**:
  - SLURMなどのジョブスケジューラ対応
  - 大規模クラスター環境での実行

- **分散実行（DistributedExecutor）**:
  - Daskベースの動的分散処理
  - クラウド環境での柔軟なスケーリング

各エグゼキュータは同一インターフェースを実装し、実行環境に応じた最適な実装を選択可能にします。

## 5. データ管理

### 5.1 永続化戦略

以下の永続化オプションをサポートします：

- **ファイルベース**:
  - CSV, JSON, YAMLでの結果保存
  - 分子ファイル（PDB, MOL2, SDF）による構造保存

- **リレーショナルデータベース**:
  - SQLAlchemyによるORM
  - 構造化クエリのサポート

- **NoSQLデータベース**:
  - MongoDBによる柔軟なドキュメント保存
  - 大規模結果セットの効率的管理

- **オブジェクトストレージ**:
  - S3互換ストレージによる大規模データ管理
  - 分散環境での共有データアクセス

### 5.2 リポジトリパターン

各ドメインオブジェクト用のリポジトリインターフェースを定義し、永続化の詳細から分離します：

```python
# 抽象リポジトリインターフェース
class ProteinRepository(ABC):
    @abstractmethod
    def save(self, protein: Protein) -> None: ...
    
    @abstractmethod
    def find_by_id(self, protein_id: str) -> Optional[Protein]: ...
    
    @abstractmethod
    def find_all(self) -> List[Protein]: ...
```

具体的な実装は各ストレージバックエンドに合わせて提供します。

## 6. パフォーマンスと分散処理

### 6.1 スケーラビリティの設計

以下の原則に基づいてスケーラビリティを確保します：

- **水平スケーリング**:
  - ステートレスな処理の分散化
  - データの分割処理（シャーディング）

- **リソース最適化**:
  - 適切なバッチサイズの動的決定
  - リソース使用効率の監視と最適化

- **キャッシング戦略**:
  - 中間結果のキャッシング
  - 頻繁に使用される分子構造のメモリキャッシュ

### 6.2 分散環境での最適化

大規模クラスター環境でのパフォーマンスを向上させるための戦略：

- **データローカリティ**:
  - 計算とデータの近接配置
  - データ転送の最小化

- **バランシング**:
  - 動的負荷分散
  - 適応的タスク割り当て

- **耐障害性**:
  - チェックポイントと再開機能
  - 障害検出と自動復旧

## 7. 開発ロードマップ

### 7.1 短期目標（0.2リリース）

- 分子表現の基本拡張（RDKit統合）
- AutoDock VinaとREStrettoの完全サポート
- 基本的なワークフロー定義フレームワーク
- シーケンシャルおよびローカル並列実行エンジン
- 最小限のCLIインターフェース

### 7.2 中期目標（0.5リリース）

- 複数のドッキングツール（gnina, diffdockなど）のサポート
- 結果解析と可視化フレームワーク
- クラスター実行エンジンの完全実装
- RESTful APIインターフェース
- カスタムスコア関数のプラグイン機構

### 7.3 長期目標（1.0リリース）

- 完全なデータ管理システム
- 高度なワークフロー定義（DAGベース）
- AIによるパラメータ最適化
- Webインターフェースと可視化ダッシュボード
- 分子動力学シミュレーションとの統合

## 8. 貢献ガイドライン

### 8.1 コード規約

- **PEP 8**: Pythonコードスタイルガイドラインの遵守
- **型アノテーション**: すべての関数とメソッドに型アノテーションを付与
- **ドキュメント**: Googleスタイルのdocstringによる文書化
- **テスト**: すべての新機能に対するテストケースの実装

### 8.2 開発プロセス

- **イシュートラッキング**: すべての新機能と修正はイシューとして登録
- **ブランチ戦略**: feature/, bugfix/, refactor/ プレフィックスによるブランチ命名
- **コードレビュー**: すべてのPRに対する少なくとも1人のレビュアー
- **CI/CD**: 自動テスト、リンター、型チェックによる品質保証

### 8.3 拡張ポイント

新機能を追加するための主要な拡張ポイント：

- **新しいドッキングツール**: `DockingToolABC`の実装
- **新しいエグゼキュータ**: `ExecutorABC`の実装
- **カスタムスコア関数**: `ScoringFunctionABC`の実装
- **データリポジトリ**: `RepositoryABC`の実装

## 9. ライセンスと謝辞

- **ライセンス**: MIT
- **引用**: 本ソフトウェアを研究で使用する場合の引用形式
- **謝辞**: 貢献者および関連プロジェクトへの謝辞

## 10. コミュニティと問い合わせ

- **リポジトリ**: GitHub リポジトリURL
- **ディスカッション**: GitHub Discussions / Slack チャンネル
- **バグ報告**: Issue Tracker
- **連絡先**: メールアドレス