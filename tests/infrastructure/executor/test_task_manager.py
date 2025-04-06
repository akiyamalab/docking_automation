from unittest.mock import MagicMock

import pytest

from docking_automation.infrastructure.executor import ExecutorABC, Task, TaskManager


class TestTaskManager:
    """TaskManagerクラスのテスト"""

    @pytest.fixture
    def sample_function(self):
        """テスト用の関数"""

        def add(a, b):
            return a + b

        return add

    @pytest.fixture
    def sample_tasks(self, sample_function):
        """テスト用のTaskインスタンスのリストを作成する"""
        return [
            Task(id="task1", function=sample_function, args={"a": 1, "b": 2}),
            Task(id="task2", function=sample_function, args={"a": 3, "b": 4}),
            Task(id="task3", function=sample_function, args={"a": 5, "b": 6}),
        ]

    @pytest.fixture
    def mock_executor(self):
        """テスト用のモックExecutorを作成する"""
        executor = MagicMock(spec=ExecutorABC)
        executor.execute.side_effect = lambda task: task.function(**task.args)
        return executor

    @pytest.mark.skip(reason="未実装のテスト")
    def test_initialization(self):
        """初期化のテスト"""
        task_manager = TaskManager()
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_add_task(self, sample_tasks):
        """タスク追加のテスト"""
        task_manager = TaskManager()
        for task in sample_tasks:
            task_manager.add_task(task)

        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_execute_all(self, sample_tasks, mock_executor):
        """すべてのタスク実行のテスト"""
        task_manager = TaskManager()
        for task in sample_tasks:
            task_manager.add_task(task)

        results = task_manager.execute_all(mock_executor)

        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_get_status(self, sample_tasks):
        """タスク状態取得のテスト"""
        task_manager = TaskManager()
        for task in sample_tasks:
            task_manager.add_task(task)

        status = task_manager.get_status()

        pass
