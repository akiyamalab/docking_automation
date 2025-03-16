from typing import Dict, Callable, Any


class Task:
    """
    実行するタスクを表すクラス。
    
    関数とその引数を保持し、実行時に関数を呼び出す。
    """
    
    def __init__(self, id: str, function: Callable, args: Dict[str, Any]):
        """
        Taskオブジェクトを初期化する。
        
        Args:
            id: タスクのID
            function: 実行する関数
            args: 関数に渡す引数の辞書
        """
        self.id = id
        self.function = function
        self.args = args
    
    def execute(self) -> Any:
        """
        タスクを実行する。
        
        Returns:
            関数の実行結果
        """
        raise NotImplementedError()