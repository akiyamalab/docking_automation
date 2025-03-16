from abc import ABC, abstractmethod
from typing import Any

from .task import Task

# インフラ
class ExecutorABC(ABC):
    """
    タスクを実行するための抽象基底クラス。
    """
    
    @abstractmethod
    def execute(self, task: Task) -> Any:
        """
        タスクを実行する。
        
        Args:
            task: 実行するタスク
            
        Returns:
            タスクの実行結果
        """
        pass