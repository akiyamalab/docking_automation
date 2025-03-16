
from .executor import ExecutorABC


class JobSchedulerExecutorABC(ExecutorABC):
    """
    ジョブスケジューラを用いて計算を行うための抽象基底クラス。
    """
    ...