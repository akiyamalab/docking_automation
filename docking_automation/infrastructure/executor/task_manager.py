from typing import List, Dict, Any

from .task import Task
from .executor import ExecutorABC


class TaskManager:
    """
    タスクを管理するクラス。
    
    複数のタスクを保持し、エグゼキュータを使用して実行する。
    """
    
    def __init__(self):
        """
        TaskManagerオブジェクトを初期化する。
        """
        self.tasks: List[Task] = []
    
    def add_task(self, task: Task) -> None:
        """
        タスクを追加する。
        
        Args:
            task: 追加するタスク
        """
        raise NotImplementedError()
    
    def execute_all(self, executor: ExecutorABC) -> List[Any]:
        """
        すべてのタスクを実行する。
        
        Args:
            executor: タスクを実行するためのエグゼキュータ
            
        Returns:
            タスクの実行結果のリスト
        """
        raise NotImplementedError()
    
    def get_status(self) -> Dict[str, str]:
        """
        各タスクの状態を取得する。
        
        Returns:
            タスクIDをキー、状態を値とする辞書
        """
        raise NotImplementedError()