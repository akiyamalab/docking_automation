from typing import List, Dict, Any, Optional, Set, Tuple, Callable
from enum import Enum

from .task import Task, TaskStatus
from .executor import ExecutorABC


# TODO: [DDD] タスク管理の実装を強化する
# - アプリケーションサービスとして実装し、ユースケースを明確に表現する
# - 依存性注入を活用し、テスト容易性を向上させる
# - トランザクション管理を導入し、タスクの一貫性を保証する
# - タスクのスケジューリングと優先順位付けの機能を追加する
# - 並列実行と同期メカニズムを強化する

class TaskManager:
    """
    タスクを管理するクラス。
    
    複数のタスクを保持し、エグゼキュータを使用して実行する。
    
    TODO: [DDD] アプリケーションサービスとしての振る舞いを強化する
    - ユースケースを明確に表現する
    - 依存性注入を活用する
    - トランザクション管理を導入する
    - ドメインイベントの発行と処理を実装する
    """
    
    def __init__(self, executor: Optional[ExecutorABC] = None):
        """
        TaskManagerオブジェクトを初期化する。
        
        TODO: [DDD] 初期化処理を強化する
        - 依存性注入を活用する
        - 設定の読み込みと検証を追加する
        - ドメインイベントの発行機能を初期化する
        
        Args:
            executor: タスクを実行するためのエグゼキュータ（指定しない場合は実行時に指定する必要がある）
        """
        self.tasks: List[Task] = []
        self.executor = executor
        self._domain_events: Set[Any] = set()  # TODO: [DDD] ドメインイベントの型を定義する
    
    def add_task(self, task: Task) -> None:
        """
        タスクを追加する。
        
        TODO: [DDD] タスク追加処理を強化する
        - タスクの重複チェックを追加する
        - 依存関係の検証を追加する
        - ドメインイベントを発行する
        
        Args:
            task: 追加するタスク
        """
        # タスクの重複チェック
        for existing_task in self.tasks:
            if existing_task.id == task.id:
                raise ValueError(f"タスクID '{task.id}' は既に存在します")
        
        self.tasks.append(task)
        # TODO: [DDD] ドメインイベントの登録機能を実装する
        # self._register_domain_event("TaskAdded", task.id)
    
    def execute_all(self, executor: Optional[ExecutorABC] = None) -> List[Any]:
        """
        すべてのタスクを実行する。
        
        TODO: [DDD] タスク実行処理を強化する
        - 依存関係に基づいた実行順序の決定を実装する
        - 並列実行のサポートを追加する
        - エラーハンドリングを強化する
        - 実行状況のモニタリングを実装する
        
        Args:
            executor: タスクを実行するためのエグゼキュータ（コンストラクタで指定していない場合は必須）
            
        Returns:
            タスクの実行結果のリスト
        """
        # エグゼキュータの取得
        actual_executor = executor or self.executor
        if actual_executor is None:
            raise ValueError("エグゼキュータが指定されていません")
        
        # 依存関係に基づいてタスクをソート
        # TODO: [DDD] 依存関係に基づいたソートを実装する
        
        # 複数のタスクを並列実行
        try:
            # execute_manyメソッドが実装されている場合は、それを使用
            if hasattr(actual_executor, 'execute_many') and callable(getattr(actual_executor, 'execute_many')):
                results = actual_executor.execute_many(self.tasks)
                # TODO: [DDD] ドメインイベントの登録機能を実装する
                # for task in self.tasks:
                #     self._register_domain_event("TaskExecuted", task.id)
            else:
                # execute_manyが実装されていない場合は、従来通り1つずつ実行
                results = []
                for task in self.tasks:
                    try:
                        # タスクを実行
                        result = actual_executor.execute(task)
                        results.append(result)
                        # TODO: [DDD] ドメインイベントの登録機能を実装する
                        # self._register_domain_event("TaskExecuted", task.id)
                    except Exception as e:
                        # エラーハンドリング
                        # TODO: [DDD] エラーハンドリングを強化する
                        results.append(None)
                        # TODO: [DDD] ドメインイベントの登録機能を実装する
                        # self._register_domain_event("TaskFailed", task.id)
        except Exception as e:
            # 並列実行中のエラーハンドリング
            print(f"タスクの並列実行中にエラーが発生しました: {e}")
            # すべてのタスクに対してNoneを返す
            results = [None] * len(self.tasks)
            # TODO: [DDD] ドメインイベントの登録機能を実装する
            # for task in self.tasks:
            #     self._register_domain_event("TaskFailed", task.id)
        
        return results
    
    def get_status(self) -> Dict[str, TaskStatus]:
        """
        各タスクの状態を取得する。
        
        TODO: [DDD] 状態取得処理を強化する
        - タスクの状態を適切に表現する
        - 状態変化の履歴を追跡する
        - 状態に基づいたフィルタリング機能を追加する
        
        Returns:
            タスクIDをキー、状態を値とする辞書
        """
        return {task.id: task.get_status() for task in self.tasks}
    
    # TODO: [DDD] ドメインイベントの登録機能を追加する
    def _register_domain_event(self, event_type: str, task_id: str) -> None:
        """
        ドメインイベントを登録する。
        
        Args:
            event_type: イベントの種類
            task_id: イベントに関連するタスクのID
        """
        # TODO: [DDD] 適切なドメインイベントクラスを実装する
        self._domain_events.add((event_type, task_id))
    
    # TODO: [DDD] ドメインイベントの取得機能を追加する
    def get_domain_events(self) -> Set[Any]:
        """
        登録されたドメインイベントを取得する。
        
        Returns:
            ドメインイベントのセット
        """
        return self._domain_events.copy()
    
    # TODO: [DDD] ドメインイベントのクリア機能を追加する
    def clear_domain_events(self) -> None:
        """
        登録されたドメインイベントをクリアする。
        """
        self._domain_events.clear()
    
    # TODO: [DDD] タスクの取得機能を追加する
    def get_task(self, task_id: str) -> Optional[Task]:
        """
        指定されたIDのタスクを取得する。
        
        Args:
            task_id: 取得するタスクのID
            
        Returns:
            タスク、見つからなければNone
        """
        for task in self.tasks:
            if task.id == task_id:
                return task
        return None
    
    # TODO: [DDD] タスクの削除機能を追加する
    def remove_task(self, task_id: str) -> bool:
        """
        指定されたIDのタスクを削除する。
        
        Args:
            task_id: 削除するタスクのID
            
        Returns:
            削除に成功した場合はTrue、そうでなければFalse
        """
        task = self.get_task(task_id)
        if task is not None:
            self.tasks.remove(task)
            # TODO: [DDD] ドメインイベントの登録機能を実装する
            # self._register_domain_event("TaskRemoved", task_id)
            return True
        return False
    
    # TODO: [DDD] タスクのフィルタリング機能を追加する
    def filter_tasks(self, condition: Callable[[Task], bool]) -> List[Task]:
        """
        条件に合致するタスクのみを取得する。
        
        Args:
            condition: フィルタリング条件
            
        Returns:
            条件に合致するタスクのリスト
        """
        return [task for task in self.tasks if condition(task)]
    
    # TODO: [DDD] タスクの状態別取得機能を追加する
    def get_tasks_by_status(self, status: TaskStatus) -> List[Task]:
        """
        指定された状態のタスクを取得する。
        
        Args:
            status: 取得するタスクの状態
            
        Returns:
            指定された状態のタスクのリスト
        """
        return self.filter_tasks(lambda task: task.get_status() == status)
    
    # TODO: [DDD] ファクトリメソッドを追加する
    @classmethod
    def create(cls, executor: Optional[ExecutorABC] = None) -> 'TaskManager':
        """
        TaskManagerオブジェクトを作成するファクトリメソッド。
        
        Args:
            executor: タスクを実行するためのエグゼキュータ（指定しない場合は実行時に指定する必要がある）
            
        Returns:
            作成されたTaskManagerオブジェクト
        """
        return cls(executor=executor)