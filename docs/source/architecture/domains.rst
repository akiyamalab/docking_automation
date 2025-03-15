==============
ドメインモデル
==============

システムは5つの主要なドメインで構成されており、それぞれが特定の責務を担当しています。
このページでは、各ドメインの概要とそのコンポーネントについて説明します。

分子準備ドメイン
===============

分子準備ドメインは、ドッキング計算に使用する分子（化合物とタンパク質）の前処理を担当します。

.. graphviz::

   digraph molecule_domain {
      rankdir=BT;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      Molecule [label="Molecule\n(抽象)"];
      Compound [label="Compound"];
      Protein [label="Protein"];
      MoleculeStructure [label="MoleculeStructure\n<<ValueObject>>"];
      MoleculeFormat [label="MoleculeFormat\n<<ValueObject>>"];
      MoleculeProperty [label="MoleculeProperty\n<<ValueObject>>"];
      
      Compound -> Molecule [arrowhead="empty", style="dashed"];
      Protein -> Molecule [arrowhead="empty", style="dashed"];
      Molecule -> MoleculeStructure [label="has"];
      Molecule -> MoleculeFormat [label="has"];
      Molecule -> MoleculeProperty [label="has many"];
   }

主要コンポーネント：

* **Molecule（抽象）**: 全ての分子の基底クラス
* **Compound**: 化合物を表現
* **Protein**: タンパク質を表現
* **MoleculeStructure**: 分子構造を表現する値オブジェクト
* **MoleculeFormat**: 分子ファイル形式を表現する値オブジェクト
* **FormatConverterService**: 分子形式変換を担当するサービス
* **MoleculePreparationService**: 分子準備を担当するサービス

ドッキング実行ドメイン
===================

ドッキング実行ドメインは、様々なドッキングツールによる計算実行を担当します。

.. graphviz::

   digraph docking_domain {
      rankdir=BT;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      DockingTask [label="DockingTask"];
      Ligand [label="Ligand"];
      Receptor [label="Receptor"];
      DockingConfiguration [label="DockingConfiguration\n<<ValueObject>>"];
      GridBox [label="GridBox\n<<ValueObject>>"];
      DockingParameter [label="DockingParameter\n<<ValueObject>>"];
      DockingToolStrategy [label="DockingToolStrategy\n<<Interface>>"];
      
      DockingTask -> Ligand [label="has"];
      DockingTask -> Receptor [label="has"];
      DockingTask -> DockingConfiguration [label="has"];
      DockingConfiguration -> GridBox [label="has"];
      DockingConfiguration -> DockingParameter [label="has many"];
      DockingTask -> DockingToolStrategy [label="uses", style="dashed"];
   }

主要コンポーネント：

* **DockingTask**: ドッキング計算タスクを表現
* **Ligand**: ドッキング用に準備された化合物
* **Receptor**: ドッキング用に準備されたタンパク質
* **DockingConfiguration**: ドッキング設定を表現する値オブジェクト
* **GridBox**: ドッキング計算の探索空間を表現する値オブジェクト
* **DockingService**: ドッキング計算の実行を担当するサービス
* **DockingToolStrategy**: 異なるドッキングツールのインターフェースを定義

ジョブ管理ドメイン
===============

ジョブ管理ドメインは、計算ジョブの分散と監視を担当します。

.. graphviz::

   digraph job_domain {
      rankdir=BT;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      Job [label="Job"];
      Task [label="Task"];
      ComputeResource [label="ComputeResource"];
      JobStatus [label="JobStatus\n<<Enumeration>>"];
      TaskStatus [label="TaskStatus\n<<Enumeration>>"];
      ResourceRequirement [label="ResourceRequirement\n<<ValueObject>>"];
      
      Job -> Task [label="has many"];
      Job -> JobStatus [label="has"];
      Job -> ComputeResource [label="uses"];
      Task -> TaskStatus [label="has"];
      Task -> DockingTask [label="references"];
      ComputeResource -> ResourceRequirement [label="has"];
   }

主要コンポーネント：

* **Job**: 計算ジョブを表現
* **Task**: ジョブ内の個別タスクを表現
* **ComputeResource**: 計算リソースを表現
* **JobStatus**: ジョブのステータスを表現する列挙型
* **ResourceRequirement**: 必要な計算リソースを表現する値オブジェクト
* **JobSchedulerService**: ジョブスケジューリングを担当するサービス
* **TaskDistributionService**: タスク分散を担当するサービス

結果管理ドメイン
=============

結果管理ドメインは、計算結果の解析と評価を担当します。

.. graphviz::

   digraph result_domain {
      rankdir=BT;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      DockingResult [label="DockingResult"];
      Pose [label="Pose\n<<ValueObject>>"];
      Score [label="Score\n<<ValueObject>>"];
      BindingInteraction [label="BindingInteraction\n<<ValueObject>>"];
      ResultAnalysis [label="ResultAnalysis"];
      ResultMetadata [label="ResultMetadata\n<<ValueObject>>"];
      
      DockingResult -> Pose [label="has many"];
      DockingResult -> Score [label="has many"];
      DockingResult -> ResultMetadata [label="has"];
      DockingResult -> DockingTask [label="references"];
      Pose -> BindingInteraction [label="has many"];
      ResultAnalysis -> DockingResult [label="analyzes"];
   }

主要コンポーネント：

* **DockingResult**: ドッキング計算結果を表現
* **Pose**: 予測された結合構造を表現する値オブジェクト
* **Score**: ドッキングスコアを表現する値オブジェクト
* **BindingInteraction**: 結合相互作用を表現する値オブジェクト
* **ResultAnalysis**: 結果解析情報を表現
* **ResultAnalysisService**: 結果解析を担当するサービス
* **VisualizationService**: 視覚化を担当するサービス

データアクセスドメイン
===================

データアクセスドメインは、データの永続化と取得を担当します。

.. graphviz::

   digraph storage_domain {
      rankdir=BT;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      DataRepository [label="DataRepository\n<<Interface>>"];
      DataRecord [label="DataRecord"];
      DataQuery [label="DataQuery\n<<ValueObject>>"];
      StorageFormat [label="StorageFormat\n<<Enumeration>>"];
      
      DataRepository -> DataRecord [label="manages"];
      DataRepository -> DataQuery [label="accepts"];
      DataRepository -> StorageFormat [label="uses"];
   }

主要コンポーネント：

* **DataRepository**: データリポジトリのインターフェース
* **DataRecord**: データレコードを表現
* **DataQuery**: データクエリを表現する値オブジェクト
* **StorageFormat**: データ保存形式を表現する列挙型
* **DataPersistenceService**: データ永続化を担当するサービス
* **DataRetrievalService**: データ取得を担当するサービス

ドメイン間の関係
==============

各ドメインは独立して設計されていますが、ドメイン間には以下のような関係があります：

.. graphviz::

   digraph domain_relationships {
      rankdir=LR;
      node [shape=box, style=filled, fillcolor=lightblue, fontname="sans-serif"];
      
      Molecule [label="分子準備ドメイン"];
      Docking [label="ドッキング実行ドメイン"];
      Job [label="ジョブ管理ドメイン"];
      Result [label="結果管理ドメイン"];
      Storage [label="データアクセスドメイン"];
      
      Molecule -> Docking [label="準備した分子を提供"];
      Docking -> Job [label="実行するタスクを提供"];
      Job -> Result [label="実行結果を提供"];
      Result -> Storage [label="結果を保存"];
      Molecule -> Storage [label="分子データを保存"];
      Docking -> Storage [label="設定を保存"];
      Job -> Storage [label="ジョブ状態を保存"];
   }

この設計により、各ドメインは独自の責務に集中し、明確な境界を持ちながらも、ドメイン間で必要な連携を行うことができます。