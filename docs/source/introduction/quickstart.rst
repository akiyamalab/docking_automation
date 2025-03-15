============
クイックスタート
============

このガイドでは、docking_automationを使った基本的なタンパク質-化合物ドッキングの実行方法を説明します。
具体的な例を通じて、システムの基本機能を理解しましょう。

サンプルデータの準備
------------

まずはサンプルデータを準備します。docking_automationにはテスト用のサンプルデータが含まれていますが、
ここでは外部から取得する方法も含めて説明します。

**タンパク質データの取得**

PDBデータベースからタンパク質構造をダウンロードする例：

.. code-block:: bash

   # PDBからタンパク質構造をダウンロード（例：2HU4）
   wget https://files.rcsb.org/download/2HU4.pdb -O protein.pdb
   
   # または以下のPythonコードでダウンロード
   python -c "import urllib.request; urllib.request.urlretrieve('https://files.rcsb.org/download/2HU4.pdb', 'protein.pdb')"

**化合物データの準備**

一般的な化合物データベースまたはサンプルファイルを使用します：

.. code-block:: bash

   # プロジェクト付属のサンプルデータの場合
   cp /path/to/docking_automation/examples/ligands/aspirin.sdf ./ligand.sdf
   
   # または一般的な医薬品のSDF/MOL形式のファイルを取得

単一分子のドッキング計算
---------------

単一のタンパク質-化合物ペアに対するドッキング計算を実行する最も基本的な方法を示します。

**Pythonコードによる実行**

.. code-block:: python

   from docking_automation.application.workflows.simple_docking_workflow import SimpleDockingWorkflow
   
   # ワークフローの作成
   workflow = SimpleDockingWorkflow()
   
   # ドッキング計算の実行
   result = workflow.execute(
       ligand_path="ligand.sdf",
       receptor_path="protein.pdb",
       config={
           "center_x": -26.75,  # ドッキング中心のX座標
           "center_y": 5.82,    # ドッキング中心のY座標
           "center_z": -8.13,   # ドッキング中心のZ座標
           "size_x": 20.0,      # ドッキングボックスのX方向サイズ
           "size_y": 20.0,      # ドッキングボックスのY方向サイズ
           "size_z": 20.0,      # ドッキングボックスのZ方向サイズ
           "exhaustiveness": 8,  # 探索の徹底度
           "num_modes": 9,       # 出力するポーズの数
       }
   )
   
   # 結果の処理
   print(f"Best score: {result.scores[0].value}")
   
   # 結果の保存
   result.save("output/docking_result.csv")

**出力例**

.. code-block:: text

   Best score: -8.4
   
   # 保存されたdocking_result.csvの内容例:
   pose_id,ligand_id,receptor_id,score,rmsd
   P001,L001,R001,-8.4,0.0
   P002,L001,R001,-7.9,1.2
   P003,L001,R001,-7.6,1.5
   ...

**コマンドラインからの実行**

.. code-block:: bash

   # 基本的なドッキング計算の実行
   docking run --ligand ligand.sdf --receptor protein.pdb --center -26.75,5.82,-8.13 --size 20,20,20 --out output/
   
   # 詳細オプションの指定
   docking run --ligand ligand.sdf --receptor protein.pdb \
       --center -26.75,5.82,-8.13 \
       --size 20,20,20 \
       --exhaustiveness 8 \
       --num-modes 9 \
       --cpu 4 \
       --out output/

**出力例**

.. code-block:: text

   Starting docking calculation...
   Preparing ligand from ligand.sdf...
   Preparing receptor from protein.pdb...
   Running docking with AutoDock Vina...
   
   Docking completed successfully.
   Found 9 binding poses.
   Best binding score: -8.4 kcal/mol
   
   Results saved to output/docking_summary.txt
   Docking poses saved to output/docking_result.pdbqt

複数化合物のドッキング計算
----------------

複数の化合物に対するバッチドッキング計算の実行方法を示します。

**Pythonコードによる実行**

.. code-block:: python

   from docking_automation.application.workflows.batch_docking_workflow import BatchDockingWorkflow
   
   # バッチワークフローの作成
   workflow = BatchDockingWorkflow()
   
   # バッチドッキング計算の実行
   ligand_paths = [
       "ligands/compound1.sdf",
       "ligands/compound2.sdf",
       "ligands/compound3.sdf"
   ]
   results = workflow.execute(
       ligand_paths=ligand_paths,  # 複数のリガンドファイルパス
       receptor_path="protein.pdb",
       output_dir="output/",
       config={
           "center_x": -26.75,
           "center_y": 5.82,
           "center_z": -8.13,
           "size_x": 20.0,
           "size_y": 20.0,
           "size_z": 20.0,
       },
       n_jobs=4  # 並列実行数
   )
   
   # 結果のランキング
   ranked_results = results.rank_by_score()
   
   # 上位10件の出力
   for i, result in enumerate(ranked_results[:10]):
       print(f"Rank {i+1}: {result.ligand.compound.id} - Score: {result.scores[0].value}")
   
   # 結果の保存
   ranked_results.save("output/ranked_results.csv")

**コマンドラインからの実行**

.. code-block:: bash

   # 複数の化合物ファイルに対するバッチドッキング
   docking batch --ligands ligands/compound1.sdf ligands/compound2.sdf ligands/compound3.sdf \
       --receptor protein.pdb \
       --center -26.75,5.82,-8.13 \
       --size 20,20,20 \
       --jobs 4 \
       --out output/

**出力例**

.. code-block:: text

   Starting batch docking with 3 ligands...
   Processing ligand 1/3: ligands/compound1.sdf
   Processing ligand 2/3: ligands/compound2.sdf
   Processing ligand 3/3: ligands/compound3.sdf
   
   Batch docking completed successfully.
   
   Summary of results:
   - compound1: Best score = -8.4 kcal/mol
   - compound2: Best score = -7.6 kcal/mol
   - compound3: Best score = -9.2 kcal/mol
   
   Results saved to output/ranked_results.csv

結果の解析と可視化
-----------

ドッキング結果の解析と可視化の基本的な方法を示します。
特に、ドッキングスコアのヒストグラム出力など、シンプルな解析機能について紹介します。

.. code-block:: python

   from docking_automation.application.services.result_analysis_service import ResultAnalysisService
   from docking_automation.domain.result.entity.docking_result import DockingResult
   import matplotlib.pyplot as plt
   
   # 結果の読み込み
   results = DockingResult.load_batch("output/ranked_results.csv")
   
   # 解析サービスの作成
   analysis_service = ResultAnalysisService()
   
   # スコアの取得
   scores = [result.get_best_score().value for result in results]
   
   # スコアのヒストグラム作成
   plt.figure(figsize=(10, 6))
   plt.hist(scores, bins=20, alpha=0.7, color='skyblue')
   plt.xlabel('ドッキングスコア (kcal/mol)')
   plt.ylabel('化合物数')
   plt.title('ドッキングスコア分布')
   plt.grid(True, alpha=0.3)
   plt.savefig('output/score_histogram.png', dpi=300, bbox_inches='tight')
   
   # 基本的な統計情報の出力
   stats = analysis_service.calculate_score_statistics(scores)
   print(f"平均スコア: {stats['mean']:.2f} kcal/mol")
   print(f"最良スコア: {stats['min']:.2f} kcal/mol")
   print(f"標準偏差: {stats['std']:.2f}")

**出力例**

.. code-block:: text

   平均スコア: -7.85 kcal/mol
   最良スコア: -9.20 kcal/mol
   標準偏差: 0.76

**ヒストグラム出力例**

.. figure:: ../images/score_histogram_example.png
   :alt: ドッキングスコアのヒストグラム例
   :width: 600px
   
   ドッキングスコア分布のヒストグラム例

次のステップ
--------

基本的な使い方を理解したら、以下のトピックを探索してより高度な機能を活用してください：

* :doc:`../api/domain/index` - APIリファレンスで利用可能な機能の詳細を確認
* :doc:`../architecture/overview` - システムのアーキテクチャを理解
* :doc:`../tutorials/index` - より詳細なチュートリアルで具体的なユースケースを学習