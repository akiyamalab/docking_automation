.. ドッキング自動化システム documentation master file

==================================
ドッキング自動化システム ドキュメント
==================================

タンパク質-化合物ドッキング計算の自動化システム（ドメイン駆動設計アプローチ）

概要
====

このプロジェクトは、タンパク質-化合物ドッキング計算を効率的に自動化するためのシステムです。ドメイン駆動設計（DDD）の原則に基づいて設計されており、以下の特徴を持ちます：

* 複数のドッキングツール（AutoDock Vina, gnina, diffdock, REstretto等）のサポート
* 柔軟なデータ保存形式（CSV, JSON, RDB, NoSQL）
* 大規模計算のための分散処理対応（ローカル実行、SLURM、クラウド）
* 明確なドメインモデルによる拡張性と保守性の高さ

主要機能
---------

* 分子準備：化合物とタンパク質の前処理
* ドッキング実行：様々なドッキングツールによる計算
* ジョブ管理：計算ジョブの分散と監視
* 結果管理：計算結果の解析と評価
* データ管理：複数形式での永続化と取得

目次
====

.. toctree::
   :maxdepth: 2
   :caption: 入門ガイド:

   introduction/overview
   introduction/installation
   introduction/quickstart

.. toctree::
   :maxdepth: 2
   :caption: アーキテクチャ:

   architecture/overview
   architecture/domains

.. toctree::
   :maxdepth: 2
   :caption: API リファレンス:

   api/domain/index
   api/application/index
   api/infrastructure/index
   api/presentation/index

.. toctree::
   :maxdepth: 1
   :caption: 開発者向け情報:

   dev/contributing

索引と検索
=========

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
