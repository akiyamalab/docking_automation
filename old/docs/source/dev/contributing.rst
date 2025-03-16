===========
開発への貢献
===========

docking_automation プロジェクトへの貢献を検討いただき、ありがとうございます。
このガイドでは、プロジェクトに貢献するための手順と注意点について説明します。

開発環境のセットアップ
---------------

1. **リポジトリのクローン**

   .. code-block:: bash

      git clone https://github.com/yourusername/docking_automation.git
      cd docking_automation

2. **開発モードでのインストール**

   .. code-block:: bash

      # 仮想環境の作成（推奨）
      python -m venv venv
      source venv/bin/activate  # Windows: venv\Scripts\activate

      # 開発モードでのインストール
      pip install -e .
      pip install -r requirements-dev.txt

3. **依存ツールのインストール**

   各種ドッキングツールや依存ライブラリのインストールについては :doc:`../introduction/installation` を参照してください。

コーディングスタイル
--------------

プロジェクトでは以下のコーディング規約に準拠しています：

* **PEP 8**: Python公式のスタイルガイド
* **Type Hints**: 静的型ヒントを使用
* **Docstrings**: Googleスタイルのドキュメント文字列

自動フォーマッティングとリンターにより、コードスタイルの一貫性を確保しています：

.. code-block:: bash

   # コードフォーマット
   black docking_automation/

   # 型チェック
   mypy docking_automation/

   # リンター
   flake8 docking_automation/

コミットメッセージ規約
--------------

コミットメッセージは以下の形式に従ってください：

.. code-block:: text

   [種別] 簡潔な説明（50文字以内）

   詳細な説明（オプション）。必要に応じて、変更の理由や
   関連する問題などを記載してください。

   Issue: #123

**種別の例**:

* **[feat]**: 新機能の追加
* **[fix]**: バグ修正
* **[docs]**: ドキュメントの変更のみ
* **[style]**: コードの動作に影響しない変更（フォーマットなど）
* **[refactor]**: リファクタリング（機能追加やバグ修正を含まない）
* **[perf]**: パフォーマンス改善
* **[test]**: テストの追加・修正
* **[chore]**: ビルドプロセスの変更など

プルリクエストの手順
-------------

1. **フォークとブランチの作成**

   GitHubでリポジトリをフォークし、新しいブランチを作成します：

   .. code-block:: bash

      git checkout -b feature/your-feature-name

2. **変更の実装**

   機能追加やバグ修正を実装します。
   変更が完了したら、テストを追加・更新し、すべてのテストがパスすることを確認します。

3. **テストの実行**

   .. code-block:: bash

      pytest

4. **コミットとプッシュ**

   変更をコミットし、フォークしたリポジトリにプッシュします：

   .. code-block:: bash

      git add .
      git commit -m "[feat] Add new feature"
      git push origin feature/your-feature-name

5. **プルリクエストの作成**

   GitHubで本リポジトリに対してプルリクエストを作成します。
   PR内容には以下を記載してください：

   * 変更の目的と概要
   * 実装したアプローチ
   * テスト方法
   * 関連するIssue

テスト
----

プロジェクトでは以下の種類のテストを使用しています：

1. **ユニットテスト**

   個々のコンポーネントの機能をテストします：

   .. code-block:: bash

      pytest docking_automation/tests/unit/

2. **統合テスト**

   複数のコンポーネントの連携をテストします：

   .. code-block:: bash

      pytest docking_automation/tests/integration/

3. **シナリオテスト**

   エンドツーエンドのワークフローをテストします：

   .. code-block:: bash

      pytest docking_automation/tests/scenario/

新しい機能やバグ修正には、適切なテストを追加してください。

ドキュメンテーション
--------------

コードの変更には、適切なドキュメントの更新が必要です：

1. **Docstrings**

   すべての公開API（クラス、メソッド、関数）には、Googleスタイルのdocstringを記述してください：

   .. code-block:: python

      def function_name(param1: type, param2: type) -> return_type:
          """関数の簡潔な説明。

          より詳細な説明を記述します。複数行にわたっても構いません。

          Args:
              param1: 引数1の説明
              param2: 引数2の説明

          Returns:
              戻り値の説明

          Raises:
              ValueError: エラーの条件を説明
          """

2. **Sphinxドキュメント**

   主要な機能や変更に関連するSphinxドキュメントも更新してください：

   .. code-block:: bash

      cd docs
      make html

   生成されたドキュメントを確認し、問題がないことを確認してください。

リリースプロセス
-----------

リリースは以下の手順で行われます：

1. バージョン番号の更新（`pyproject.toml`）
2. CHANGELOG.mdの更新
3. リリースブランチの作成
4. CIパイプラインによるテストとビルド
5. リリースタグの作成
6. PyPIへのパッケージのアップロード

質問やサポート
----------

質問や支援が必要な場合は、以下の方法で連絡できます：

* GitHubの Issue を作成
* 開発チームのメーリングリストにメッセージを送信
* コミュニティチャットに参加

よくある質問や既知の問題については :doc:`troubleshooting` を参照してください。