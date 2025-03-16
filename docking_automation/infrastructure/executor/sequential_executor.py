
from .executor import ExecutorABC


# インフラ
class SequentialExecutor(ExecutorABC):
    """
    逐次実行を行うクラス
    """
    ...