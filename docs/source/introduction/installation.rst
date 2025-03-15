==========
インストール
==========

システム要件
--------

docking_automationを使用するためには、以下の要件を満たす必要があります：

**ソフトウェア要件**
  * **Python**: バージョン3.8以上
  * **外部ツール**: 各種ドッキングツール（下記参照）

依存ツールのインストール
---------------

docking_automationは、以下の外部ツールを使用します。各ツールのインストール方法は、それぞれの公式ドキュメントを参照してください。

**必須ツール**

* **Open Babel**: 分子ファイル形式の変換
  
  .. code-block:: bash
  
     # Ubuntu/Debian
     sudo apt-get install openbabel
     
     # Anaconda
     conda install -c conda-forge openbabel

* **AutoDock Vina**: 基本ドッキングツール
  
  .. code-block:: bash
  
     # 手動インストール
     wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz
     tar -xzf autodock_vina_1_1_2_linux_x86.tgz
     export PATH=$PATH:$PWD/autodock_vina_1_1_2_linux_x86/bin
     
     # Anaconda
     conda install -c bioconda autodock-vina

* **MGLTools**: 分子準備ツール (prepare_receptor, prepare_ligand)
  
  MGLToolsは公式サイト(https://ccsb.scripps.edu/mgltools/)からダウンロードしてインストールしてください。

**オプションツール**

* **gnina**: 機械学習ベースのドッキングツール
* **diffdock**: 拡散モデルを用いたドッキングツール
* **REstretto**: レセプターアンサンブルを用いたドッキングツール

パッケージのインストール
---------------

docking_automationは、以下の方法でインストールできます：

PyPIからのインストール（推奨）
^^^^^^^^^^^^^^^

安定版をPyPIからインストールするには：

.. code-block:: bash

   pip install docking_automation

これにより、最新の安定版が自動的にインストールされます。

特定のバージョンを指定する場合：

.. code-block:: bash

   pip install docking_automation==0.1.0

Anacondaでのインストール（推奨）
^^^^^^^^^^^^^^^^^^^

Anacondaユーザー向けのインストール方法：

.. code-block:: bash

   # 仮想環境の作成（推奨）
   conda create -n docking-env python=3.8
   conda activate docking-env
   
   # 依存パッケージのインストール
   conda install -c conda-forge rdkit openbabel
   
   # docking_automationのインストール
   pip install docking_automation

開発版としてインストール
^^^^^^^^^^^^^^^^^^^^^^^

最新の開発版を使用するには、Gitリポジトリからクローンし、開発モードでインストールします：

.. code-block:: bash

   git clone https://github.com/yourusername/docking_automation.git
   cd docking_automation
   pip install -e .

このインストール方法では、ソースコードを編集しても再インストールせずに変更が反映されます。

インストールの確認
------------

インストールが正常に完了したことを確認するには、以下のコマンドを実行します：

.. code-block:: bash

   # Python APIのインポートテスト
   python -c "from docking_automation.application.workflows.simple_docking_workflow import SimpleDockingWorkflow; print('インストール成功')"
   
   # CLIツールの確認
   docking --version
