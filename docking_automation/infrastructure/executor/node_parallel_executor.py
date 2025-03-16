
from .executor import ExecutorABC


# インフラ
class NodeParallelExecutor(ExecutorABC):
    """
    ノード内並列実行を行うクラス
    """
    ...