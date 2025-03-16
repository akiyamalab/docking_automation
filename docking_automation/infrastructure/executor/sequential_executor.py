from typing import Any

from .executor import ExecutorABC
from .task import Task


# インフラ
class SequentialExecutor(ExecutorABC):
    """
    タスクを逐次実行するクラス。
    """
    
    def execute(self, task: Task) -> Any:
        """
        タスクを実行する。
        
        Args:
            task: 実行するタスク
            
        Returns:
            タスクの実行結果
        """
        raise NotImplementedError()