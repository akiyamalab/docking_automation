================
プレゼンテーション層API
================

プレゼンテーション層は、ユーザーとシステムの対話を担当する層です。
コマンドラインインターフェース（CLI）、プログラミングAPI、Webインターフェースなどを提供します。

役割と責務
--------

プレゼンテーション層の主な責務は以下の通りです：

* ユーザー入力の受け取りとバリデーション
* アプリケーション層の機能呼び出し
* 結果の表示とフォーマット
* エラーハンドリングとユーザーへのフィードバック

コマンドラインインターフェース (CLI)
---------------------------

.. automodule:: docking_automation.presentation.cli
   :members:
   :undoc-members:
   :show-inheritance:

docking_cli
^^^^^^^^^

メインのCLIツールの実装です。

.. automodule:: docking_automation.presentation.cli.docking_cli
   :members:
   :undoc-members:
   :show-inheritance:

CLIの基本的な使い方：

.. code-block:: bash

   # ヘルプの表示
   docking --help
   
   # 単一ドッキング計算の実行
   docking run --ligand ligand.sdf --receptor protein.pdb --center -26.75,5.82,-8.13 --size 20,20,20
   
   # バッチドッキング計算の実行
   docking batch --ligands-dir ligands/ --receptor protein.pdb --out results/
   
   # 結果の表示
   docking result show --file results/docking_results.csv

プログラミングAPI
------------

.. automodule:: docking_automation.presentation.api
   :members:
   :undoc-members:
   :show-inheritance:

RESTful API
^^^^^^^^^^

RESTful APIの実装（オプション機能）です。

.. code-block:: python

   # APIサーバーの起動例
   from docking_automation.presentation.api.rest_api import create_app
   
   app = create_app()
   app.run(host='0.0.0.0', port=5000)

API仕様：

.. list-table::
   :header-rows: 1
   :widths: 15 15 70

   * - エンドポイント
     - メソッド
     - 説明
   * - /api/v1/docking
     - POST
     - ドッキング計算を実行
   * - /api/v1/docking/{id}
     - GET
     - ドッキング結果を取得
   * - /api/v1/jobs
     - GET
     - ジョブ一覧を取得
   * - /api/v1/jobs/{id}
     - GET
     - ジョブの詳細を取得

プログラムによるAPI利用例：

.. code-block:: python

   import requests
   import json
   
   # ドッキング計算のリクエスト
   response = requests.post(
       'http://localhost:5000/api/v1/docking',
       json={
           'ligand': 'base64エンコードされたリガンドファイル',
           'receptor': 'base64エンコードされたレセプターファイル',
           'center_x': -26.75,
           'center_y': 5.82,
           'center_z': -8.13,
           'size_x': 20.0,
           'size_y': 20.0,
           'size_z': 20.0
       }
   )
   
   # 結果の処理
   result = json.loads(response.text)
   job_id = result['job_id']
   
   # ジョブ状態の確認
   job_response = requests.get(f'http://localhost:5000/api/v1/jobs/{job_id}')
   job_status = json.loads(job_response.text)
   
   print(f"Job status: {job_status['status']}")

Webインターフェース
--------------

.. automodule:: docking_automation.presentation.web
   :members:
   :undoc-members:
   :show-inheritance:

Webベースの可視化インターフェース（オプション機能）です。

機能一覧：

* ドッキング計算の設定と実行
* 結果の3D可視化
* 化合物ライブラリの管理
* 計算ジョブのモニタリング

インターフェース連携
--------------

プレゼンテーション層では、複数のインターフェースを統一的に扱うためのファクトリーパターンを使用しています。

.. automodule:: docking_automation.presentation.cli.mock_factory
   :members:
   :undoc-members:
   :show-inheritance:

エラーハンドリング
------------

プレゼンテーション層では、以下の種類のエラーハンドリングを行います：

* 入力エラー：ユーザー入力の検証
* アプリケーションエラー：アプリケーション層からのエラー
* システムエラー：外部システムやリソースに関するエラー

エラーメッセージの例：

.. code-block:: text

   ERROR: Invalid input file format. Expected SDF or MOL format for ligand.
   ERROR: Receptor file not found: receptor.pdb
   ERROR: Docking calculation failed: Unable to connect to Vina executable.
   
   # デバッグモードでの詳細なエラー情報
   DEBUG: Traceback (most recent call last):
     File "docking_automation/presentation/cli/docking_cli.py", line 123, in run_docking
       result = workflow.execute(ligand_path, receptor_path, config)
     ...