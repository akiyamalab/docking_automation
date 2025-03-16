
from .job_scheduler_executor import JobSchedulerExecutorABC


class SlurmExecutor(JobSchedulerExecutorABC):
    """
    Slurmを用いてノード間並列計算を行うためのクラス
    """
    ...