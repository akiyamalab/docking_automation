================
アプリケーション層API
================

アプリケーション層は、ドメイン層のビジネスロジックを調整し、ユースケースを実装する層です。
この層は、プレゼンテーション層とドメイン層の間の橋渡しの役割を果たします。

役割と責務
--------

アプリケーション層の主な責務は以下の通りです：

* ユースケースの実装
* ユーザー操作の調整
* ドメインサービスの連携
* データの変換と提供
* トランザクション管理

ワークフロー
---------

.. automodule:: docking_automation.application.workflows
   :members:
   :undoc-members:
   :show-inheritance:

.. autoclass:: docking_automation.application.workflows.simple_docking_workflow.SimpleDockingWorkflow
   :members:
   :undoc-members:
   :show-inheritance:
   :special-members: __init__

   .. automethod:: execute

バッチドッキングのワークフロー
^^^^^^^^^^^^^^^^^^^^^

バッチドッキングワークフローは、複数のリガンドに対して効率的にドッキング計算を実行するためのワークフローです。

.. code-block:: python

   # バッチドッキングの例
   from docking_automation.application.workflows.batch_docking_workflow import BatchDockingWorkflow
   
   workflow = BatchDockingWorkflow()
   results = workflow.execute(
       ligands_dir="ligands/",
       receptor_path="receptor.pdb",
       output_dir="results/",
       n_jobs=4
   )

アプリケーションサービス
-----------------

.. automodule:: docking_automation.application.services
   :members:
   :undoc-members:
   :show-inheritance:

主なサービスには以下があります：

* **DockingApplicationService**: ドッキング実行の調整を担当
* **MoleculePreparationApplicationService**: 分子準備の調整を担当
* **ResultAnalysisApplicationService**: 結果解析の調整を担当
* **JobManagementApplicationService**: ジョブ管理の調整を担当

データ転送オブジェクト (DTO)
----------------------

データ転送オブジェクト（DTO）は、レイヤー間でデータを交換するために使用されるオブジェクトです。

.. automodule:: docking_automation.application.dto
   :members:
   :undoc-members:
   :show-inheritance:

主なDTOには以下があります：

* **DockingRequestDTO**: ドッキング計算リクエストを表現
* **DockingResultDTO**: ドッキング計算結果を表現
* **MoleculeDTO**: 分子情報を表現
* **JobDTO**: ジョブ情報を表現

マッパー
------

マッパーは、DTOとドメインオブジェクト間の変換を担当します。

.. automodule:: docking_automation.application.mapper
   :members:
   :undoc-members:
   :show-inheritance:

インターフェース
-----------

アプリケーション層は、プレゼンテーション層に対して以下のインターフェースを提供します：

.. automodule:: docking_automation.application.interfaces
   :members:
   :undoc-members:
   :show-inheritance:

使用例
-----

基本的なドッキング計算：

.. code-block:: python

   from docking_automation.application.services.docking_application_service import DockingApplicationService
   from docking_automation.application.dto.docking_request_dto import DockingRequestDTO
   
   # アプリケーションサービスの取得
   service = DockingApplicationService()
   
   # リクエストの作成
   request = DockingRequestDTO(
       ligand_path="ligand.sdf",
       receptor_path="receptor.pdb",
       center_x=-26.75,
       center_y=5.82,
       center_z=-8.13,
       size_x=20.0,
       size_y=20.0,
       size_z=20.0,
       exhaustiveness=8
   )
   
   # ドッキング計算の実行
   result = service.execute_docking(request)
   
   # 結果の処理
   print(f"Best score: {result.scores[0]}")
   result.save_to_file("output.csv")