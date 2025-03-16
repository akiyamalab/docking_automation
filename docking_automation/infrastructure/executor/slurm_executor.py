
from .job_scheduler_executor import JobSchedulerExecutorABC


# インフラ
class SlurmExecutor(JobSchedulerExecutorABC):
    """
    Slurmを用いてノード間並列計算を行うためのクラス
    """
    ...