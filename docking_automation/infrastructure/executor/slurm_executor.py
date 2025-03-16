from typing import Any, Dict

from .job_scheduler_executor import JobSchedulerExecutorABC
from .task import Task


# インフラ
class SlurmExecutor(JobSchedulerExecutorABC):
    """
    SLURMを用いてノード間並列計算を行うためのクラス。
    """
    
    def __init__(self, slurm_config: Dict[str, Any]):
        """
        SlurmExecutorオブジェクトを初期化する。
        
        Args:
            slurm_config: SLURMの設定
        """
        self.slurm_config = slurm_config
    
    def execute(self, task: Task) -> Any:
        """
        タスクを実行する。
        
        Args:
            task: 実行するタスク
            
        Returns:
            タスクの実行結果
        """
        raise NotImplementedError()
    
    def submit_job(self, task: Task) -> str:
        """
        ジョブを投入する。
        
        Args:
            task: 投入するタスク
            
        Returns:
            ジョブID
        """
        raise NotImplementedError()
    
    def check_status(self, job_id: str) -> str:
        """
        ジョブの状態を確認する。
        
        Args:
            job_id: ジョブID
            
        Returns:
            ジョブの状態
        """
        raise NotImplementedError()