
from .job_scheduler_executor import JobSchedulerExecutorABC


class TSUBAMEUGEExecutor(JobSchedulerExecutorABC):
    """
    TSUBAME4.0 の Univa Grid Engine を用いてノード間並列計算を行うためのクラス
    """
    ...