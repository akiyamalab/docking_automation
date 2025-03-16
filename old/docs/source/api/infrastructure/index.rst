====================
インフラストラクチャ層API
====================

インフラストラクチャ層は、外部システムやフレームワークとの連携を担当する層です。
データベースアクセス、外部ツール連携、ファイル操作などの具体的な実装を提供します。

役割と責務
--------

インフラストラクチャ層の主な責務は以下の通りです：

* ドメイン層で定義されたリポジトリの実装
* 外部ツールとの連携
* データの永続化
* ファイル形式の変換
* 計算リソースの管理

外部ツール連携
---------

.. automodule:: docking_automation.infrastructure.tools
   :members:
   :undoc-members:
   :show-inheritance:

Vinaツール連携
^^^^^^^^^^

AutoDock Vinaとの連携を担当するコンポーネントです。

.. automodule:: docking_automation.infrastructure.tools.vina
   :members:
   :undoc-members:
   :show-inheritance:

.. code-block:: python

   # Vinaアダプターの使用例
   from docking_automation.infrastructure.tools.vina.vina_adapter import VinaAdapter
   
   vina_adapter = VinaAdapter(vina_executable_path="/usr/local/bin/vina")
   result = vina_adapter.execute(docking_task)

その他のドッキングツール
^^^^^^^^^^^^^^

その他のサポートされているドッキングツールのアダプターには以下があります：

* **GninaAdapter**: Gninaツールとの連携
* **DiffDockAdapter**: DiffDockツールとの連携
* **RestrettoAdapter**: RESTRettoツールとの連携

ファイル形式変換
-----------

.. automodule:: docking_automation.infrastructure.formats
   :members:
   :undoc-members:
   :show-inheritance:

OpenBabel形式変換
^^^^^^^^^^^^^

OpenBabelを使用した分子ファイル形式の変換を担当します。

.. automodule:: docking_automation.infrastructure.formats.openbabel_format_converter
   :members:
   :undoc-members:
   :show-inheritance:

.. code-block:: python

   # OpenBabel形式変換の使用例
   from docking_automation.infrastructure.formats.openbabel_format_converter import OpenBabelFormatConverter
   
   converter = OpenBabelFormatConverter()
   converted_file = converter.convert("molecule.sdf", "pdbqt")

データ永続化
---------

.. automodule:: docking_automation.infrastructure.persistence
   :members:
   :undoc-members:
   :show-inheritance:

CSVデータ永続化
^^^^^^^^^^^

CSV形式でのデータ永続化を担当します。

.. automodule:: docking_automation.infrastructure.persistence.csv
   :members:
   :undoc-members:
   :show-inheritance:

JSONデータ永続化
^^^^^^^^^^^

JSON形式でのデータ永続化を担当します。

.. automodule:: docking_automation.infrastructure.persistence.json
   :members:
   :undoc-members:
   :show-inheritance:

計算リソース管理
-----------

.. automodule:: docking_automation.infrastructure.compute
   :members:
   :undoc-members:
   :show-inheritance:

ローカル計算管理
^^^^^^^^^^^

ローカルマシンでの計算実行を担当します。

.. automodule:: docking_automation.infrastructure.compute.local
   :members:
   :undoc-members:
   :show-inheritance:

SLURM計算管理
^^^^^^^^^^

SLURMクラスターでの計算実行を担当します。

.. automodule:: docking_automation.infrastructure.compute.slurm
   :members:
   :undoc-members:
   :show-inheritance:

.. code-block:: python

   # SLURMジョブ管理の使用例
   from docking_automation.infrastructure.compute.slurm.slurm_job_manager import SlurmJobManager
   
   job_manager = SlurmJobManager()
   job_id = job_manager.submit_job(
       script_path="docking_job.sh",
       cpu_cores=4,
       memory_gb=8,
       time_limit_minutes=120
   )
   status = job_manager.get_job_status(job_id)

クラウド計算管理
^^^^^^^^^^

クラウド環境での計算実行を担当します。

.. automodule:: docking_automation.infrastructure.compute.cloud
   :members:
   :undoc-members:
   :show-inheritance:

設定管理
------

.. automodule:: docking_automation.infrastructure.config
   :members:
   :undoc-members:
   :show-inheritance:

環境設定
^^^^^^

環境変数やシステム設定の管理を担当します。

.. automodule:: docking_automation.infrastructure.config.environment
   :members:
   :undoc-members:
   :show-inheritance:

リポジトリ実装
----------

ドメイン層で定義されたリポジトリインターフェースの具体的な実装を提供します。

.. code-block:: python

   # リポジトリ実装の使用例
   from docking_automation.infrastructure.persistence.csv.csv_docking_result_repository import CsvDockingResultRepository
   
   repository = CsvDockingResultRepository(base_dir="./results")
   repository.save(docking_result)
   result = repository.find_by_id("result-123")