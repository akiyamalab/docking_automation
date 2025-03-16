import pytest
from unittest.mock import patch, MagicMock
from docking_automation.infrastructure.executor.task import Task
from docking_automation.infrastructure.executor.slurm_executor import SlurmExecutor


class TestSlurmExecutor:
    """SlurmExecutorクラスのテスト"""
    
    @pytest.fixture
    def executor(self):
        """テスト用のSlurmExecutorインスタンスを作成する"""
        return SlurmExecutor(slurm_config={
            "partition": "test",
            "time": "1:00:00",
            "mem": "4G",
            "cpus-per-task": 1
        })
    
    @pytest.fixture
    def sample_function(self):
        """テスト用の関数"""
        def add(a, b):
            return a + b
        return add
    
    @pytest.fixture
    def sample_task(self, sample_function):
        """テスト用のTaskインスタンスを作成する"""
        return Task(
            id="task1",
            function=sample_function,
            args={"a": 1, "b": 2}
        )
    
    @patch("subprocess.run")
    def test_submit_job(self, mock_run, executor, sample_task):
        """ジョブ投入のテスト"""
        # sbatchコマンドの実行結果をモック
        mock_process = MagicMock()
        mock_process.stdout = "Submitted batch job 12345"
        mock_run.return_value = mock_process
        
        job_id = executor.submit_job(sample_task)
        
        raise NotImplementedError()
    
    @patch("subprocess.run")
    def test_check_status(self, mock_run, executor):
        """ジョブ状態確認のテスト"""
        # sacctコマンドの実行結果をモック
        mock_process = MagicMock()
        mock_process.stdout = "JobID|State\n12345|RUNNING"
        mock_run.return_value = mock_process
        
        status = executor.check_status("12345")
        
        raise NotImplementedError()
    
    @patch("subprocess.run")
    def test_execute(self, mock_run, executor, sample_task):
        """タスク実行のテスト"""
        # sbatchコマンドの実行結果をモック
        mock_process_submit = MagicMock()
        mock_process_submit.stdout = "Submitted batch job 12345"
        
        # sacctコマンドの実行結果をモック
        mock_process_status = MagicMock()
        mock_process_status.stdout = "JobID|State\n12345|COMPLETED"
        
        # 結果ファイルの読み込み結果をモック
        mock_process_result = MagicMock()
        mock_process_result.stdout = "3"  # 1 + 2 = 3
        
        # 複数回のsubprocess.runの呼び出しに対して異なる戻り値を設定
        mock_run.side_effect = [mock_process_submit, mock_process_status, mock_process_result]
        
        result = executor.execute(sample_task)
        
        raise NotImplementedError()