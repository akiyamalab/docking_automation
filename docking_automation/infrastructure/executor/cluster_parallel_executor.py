
from .executor import ExecutorABC


class ClusterParallelExecutor(ExecutorABC):
    """
    ノード間並列実行を行うクラス
    """
    ...