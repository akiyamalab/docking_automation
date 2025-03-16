import pytest
from docking_automation.infrastructure.executor.task import Task
from docking_automation.infrastructure.executor.sequential_executor import SequentialExecutor


class TestSequentialExecutor:
    """SequentialExecutorクラスのテスト"""
    
    @pytest.fixture
    def executor(self):
        """テスト用のSequentialExecutorインスタンスを作成する"""
        return SequentialExecutor()
    
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
    
    def test_execute(self, executor, sample_task):
        """タスク実行のテスト"""
        result = executor.execute(sample_task)
        raise NotImplementedError()
    
    def test_execute_with_exception(self, executor):
        """例外が発生するタスク実行のテスト"""
        def failing_function():
            raise ValueError("Test error")
        
        task = Task(
            id="failing_task",
            function=failing_function,
            args={}
        )
        
        raise NotImplementedError()