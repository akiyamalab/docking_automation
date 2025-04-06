from abc import ABC, abstractmethod
from typing import Any, Dict, List

from .task import Task

# TODO: [DDD] エグゼキュータの実装を強化する
# - ポートとアダプターパターンを適用し、ドメインロジックとインフラストラクチャを分離する
# - 依存性注入を活用し、テスト容易性を向上させる
# - 実行戦略パターンを導入し、異なる実行方法を柔軟に切り替えられるようにする
# - エラー処理と再試行メカニズムを強化する
# - モニタリングと診断機能を追加する


# インフラ
class ExecutorABC(ABC):
    """
    タスクを実行するための抽象基底クラス。

    TODO: [DDD] インターフェースを強化する
    - 実行戦略パターンを導入する
    - エラー処理と再試行メカニズムを定義する
    - モニタリングと診断機能を定義する
    - 並列実行のサポートを追加する
    """

    @abstractmethod
    def execute(self, task: Task) -> Any:
        """
        タスクを実行する。

        TODO: [DDD] 実行処理を強化する
        - タスクの状態管理を実装する
        - エラー処理と再試行メカニズムを追加する
        - 実行結果の標準化を行う
        - ドメインイベントの発行をサポートする

        Args:
            task: 実行するタスク

        Returns:
            タスクの実行結果
        """
        pass

    # TODO: [DDD] 複数タスクの実行機能を追加する
    @abstractmethod
    def execute_many(self, tasks: List[Task]) -> List[Any]:
        """
        複数のタスクを実行する。

        Args:
            tasks: 実行するタスクのリスト

        Returns:
            タスクの実行結果のリスト
        """
        pass

    # TODO: [DDD] 非同期実行機能を追加する
    @abstractmethod
    def execute_async(self, task: Task) -> str:
        """
        タスクを非同期で実行する。

        Args:
            task: 実行するタスク

        Returns:
            非同期実行のジョブID
        """
        pass

    # TODO: [DDD] 実行状態の取得機能を追加する
    @abstractmethod
    def get_execution_status(self, job_id: str) -> Dict[str, Any]:
        """
        非同期実行の状態を取得する。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            実行状態を表す辞書
        """
        pass

    # TODO: [DDD] 実行結果の取得機能を追加する
    @abstractmethod
    def get_execution_result(self, job_id: str) -> Any:
        """
        非同期実行の結果を取得する。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            タスクの実行結果
        """
        pass

    # TODO: [DDD] 実行のキャンセル機能を追加する
    @abstractmethod
    def cancel_execution(self, job_id: str) -> bool:
        """
        非同期実行をキャンセルする。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            キャンセルに成功した場合はTrue、そうでなければFalse
        """
        pass

    # TODO: [DDD] ファクトリメソッドを追加する
    @classmethod
    def create(cls, executor_type: str, **kwargs) -> "ExecutorABC":
        """
        エグゼキュータを作成するファクトリメソッド。

        Args:
            executor_type: エグゼキュータの種類
            **kwargs: エグゼキュータの設定

        Returns:
            作成されたエグゼキュータ
        """
        # TODO: [DDD] エグゼキュータの種類に応じたインスタンス生成を実装する
        raise NotImplementedError("ファクトリメソッドを実装する必要があります")
