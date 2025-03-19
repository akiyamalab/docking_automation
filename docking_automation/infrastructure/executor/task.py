from typing import Dict, Callable, Any, Optional, List, Set, Union
from enum import Enum, auto


# TODO: [DDD] タスクの実装を強化する
# - 値オブジェクトとしての特性を強化し、不変性を確保する
# - タスクの状態管理を追加し、ライフサイクルを明確にする
# - タスク間の依存関係を表現できるようにする
# - タスクの実行結果を表す値オブジェクトを導入する
# - ドメインイベントを発行し、タスクの状態変化を通知できるようにする

# TODO: [DDD] タスクの状態を表す列挙型を追加する
class TaskStatus(Enum):
    """
    タスクの状態を表す列挙型。
    """
    PENDING = auto()
    RUNNING = auto()
    COMPLETED = auto()
    FAILED = auto()
    CANCELED = auto()


class Task:
    """
    実行するタスクを表すクラス。
    
    関数とその引数を保持し、実行時に関数を呼び出す。
    
    TODO: [DDD] 値オブジェクトとしての振る舞いを強化する
    - 不変性を確保する（dataclass(frozen=True)の使用を検討）
    - 等価性比較を実装する
    - 状態変化を追跡する
    - ライフサイクルを明確に定義する
    """
    
    def __init__(
        self,
        id: str,
        function: Callable,
        args: Dict[str, Any],
        dependencies: Optional[List[str]] = None,
        status: TaskStatus = TaskStatus.PENDING
    ):
        """
        Taskオブジェクトを初期化する。
        
        TODO: [DDD] 初期化処理を強化する
        - IDの一意性を保証する仕組みを導入する
        - 依存関係の検証を追加する
        - 不変条件のバリデーションを追加する
        
        Args:
            id: タスクのID
            function: 実行する関数
            args: 関数に渡す引数の辞書
            dependencies: 依存するタスクのIDのリスト
            status: タスクの初期状態
        """
        self.id = id
        self.function = function
        self.args = args
        self.dependencies = dependencies or []
        self.status = status
        self._result: Any = None
        self._error: Optional[Exception] = None
        self._domain_events: Set[Any] = set()  # TODO: [DDD] ドメインイベントの型を定義する
    
    def execute(self) -> Any:
        """
        タスクを実行する。
        
        TODO: [DDD] 実行処理を強化する
        - 状態遷移を管理する
        - エラーハンドリングを強化する
        - 実行結果を適切に保持する
        - ドメインイベントを発行する
        
        Returns:
            関数の実行結果
        """
        try:
            self.status = TaskStatus.RUNNING
            self._register_domain_event("TaskStarted")
            
            self._result = self.function(**self.args)
            self.status = TaskStatus.COMPLETED
            self._register_domain_event("TaskCompleted")
            
            return self._result
        except Exception as e:
            self.status = TaskStatus.FAILED
            self._error = e
            self._register_domain_event("TaskFailed")
            raise
    
    # TODO: [DDD] ドメインイベントの登録機能を追加する
    def _register_domain_event(self, event_type: str) -> None:
        """
        ドメインイベントを登録する。
        
        Args:
            event_type: イベントの種類
        """
        # TODO: [DDD] 適切なドメインイベントクラスを実装する
        self._domain_events.add((event_type, self.id))
    
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
    
    # TODO: [DDD] タスクの状態取得機能を追加する
    def get_status(self) -> TaskStatus:
        """
        タスクの状態を取得する。
        
        Returns:
            タスクの状態
        """
        return self.status
    
    # TODO: [DDD] タスクの結果取得機能を追加する
    def get_result(self) -> Any:
        """
        タスクの実行結果を取得する。
        
        Returns:
            タスクの実行結果（実行前または失敗した場合はNone）
        """
        return self._result
    
    # TODO: [DDD] タスクのエラー取得機能を追加する
    def get_error(self) -> Optional[Exception]:
        """
        タスクの実行エラーを取得する。
        
        Returns:
            タスクの実行エラー（エラーが発生していない場合はNone）
        """
        return self._error
    
    # TODO: [DDD] タスクのキャンセル機能を追加する
    def cancel(self) -> bool:
        """
        タスクをキャンセルする。
        
        Returns:
            キャンセルに成功した場合はTrue、そうでなければFalse
        """
        if self.status in [TaskStatus.PENDING, TaskStatus.RUNNING]:
            self.status = TaskStatus.CANCELED
            self._register_domain_event("TaskCanceled")
            return True
        return False
    
    # TODO: [DDD] エンティティの等価性比較を実装する
    def __eq__(self, other: object) -> bool:
        """
        等価性比較を行う。タスクの場合はIDのみで比較する。
        
        Args:
            other: 比較対象のオブジェクト
            
        Returns:
            等価であればTrue、そうでなければFalse
        """
        if not isinstance(other, Task):
            return False
        return self.id == other.id
    
    # TODO: [DDD] エンティティのハッシュ値計算を実装する
    def __hash__(self) -> int:
        """
        ハッシュ値を計算する。タスクの場合はIDのみを使用する。
        
        Returns:
            ハッシュ値
        """
        return hash(self.id)
    
    # TODO: [DDD] ファクトリメソッドを追加する
    @classmethod
    def create(
        cls,
        function: Callable,
        args: Dict[str, Any],
        dependencies: Optional[List[str]] = None,
        id: Optional[str] = None
    ) -> 'Task':
        """
        Taskオブジェクトを作成するファクトリメソッド。
        
        Args:
            function: 実行する関数
            args: 関数に渡す引数の辞書
            dependencies: 依存するタスクのIDのリスト
            id: タスクのID（指定しない場合は自動生成）
            
        Returns:
            作成されたTaskオブジェクト
        """
        import uuid
        task_id = id or f"task_{uuid.uuid4()}"
        return cls(
            id=task_id,
            function=function,
            args=args,
            dependencies=dependencies
        )