
from .executor import ExecutorABC

# インフラ
class ClusterParallelExecutor(ExecutorABC):
    """
    ノード間並列実行を行うクラス
    """
    ...