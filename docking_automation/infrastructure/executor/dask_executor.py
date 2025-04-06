"""
Daskを使った並列処理Executorモジュール。

このモジュールは、Daskを使って複数のタスクを並列実行するためのExecutorを提供します。
ローカル環境とクラスタ環境（Slurm、PBS等）の両方に対応しています。
"""

import uuid
from typing import Any, Dict, List, Optional, Tuple

from dask.distributed import Client, LocalCluster
from dask_jobqueue.pbs import PBSCluster
from dask_jobqueue.slurm import SLURMCluster
from tqdm import tqdm

from .executor import ExecutorABC
from .task import Task, TaskStatus


class DaskExecutor(ExecutorABC):
    """
    Daskを使って並列処理を行うExecutor。

    複数のDockingToolABC.run_docking()メソッド呼び出しを並列に実行します。
    ローカル環境とクラスタ環境（Slurm、PBS等）の両方に対応しています。

    Attributes:
        scheduler_type: スケジューラの種類 ("local", "slurm", "pbs", "sge" など)
        n_workers: ローカル環境での並列ワーカー数
        scheduler_kwargs: スケジューラ固有の追加パラメータ
    """

    def __init__(self, scheduler_type: str = "local", n_workers: Optional[int] = None, **scheduler_kwargs):
        """
        DaskExecutorオブジェクトを初期化します。

        Args:
            scheduler_type: スケジューラの種類 ("local", "slurm", "pbs", "sge" など)
            n_workers: ローカル環境での並列ワーカー数（Noneの場合は自動設定）
            **scheduler_kwargs: スケジューラ固有の追加パラメータ
        """
        self.scheduler_type = scheduler_type
        self.n_workers = n_workers
        self.scheduler_kwargs = scheduler_kwargs
        self._client = None
        self._futures: Dict[str, Tuple[Any, Task]] = {}

    def _setup_client(self):
        """
        環境に応じたDaskクライアントをセットアップします。

        Returns:
            設定されたDaskクライアント
        """

        if self.scheduler_type == "local":
            # ローカル環境

            cluster = LocalCluster(n_workers=self.n_workers)
            client = Client(cluster)
            print(f"ローカル環境で実行します（ワーカー数: {client.ncores}）")
            return client
        elif self.scheduler_type == "slurm":

            cluster = SLURMCluster(**self.scheduler_kwargs)
            cluster.scale(jobs=self.scheduler_kwargs.get("jobs", 10))
            client = Client(cluster)
            print(f"Slurm環境で実行します（ジョブ数: {self.scheduler_kwargs.get('jobs', 10)}）")
            return client
        elif self.scheduler_type == "pbs":
            # PBS環境
            cluster = PBSCluster(**self.scheduler_kwargs)
            cluster.scale(jobs=self.scheduler_kwargs.get("jobs", 10))
            client = Client(cluster)
            print(f"PBS環境で実行します（ジョブ数: {self.scheduler_kwargs.get('jobs', 10)}）")
            return client
        else:
            raise ValueError(f"未対応のスケジューラタイプです: {self.scheduler_type}")

    def execute(self, task: Task) -> Any:
        """
        タスクを実行します。

        Args:
            task: 実行するタスク

        Returns:
            タスクの実行結果
        """
        # 単一タスクの場合は直接実行
        return task.execute()

    def execute_many(self, tasks: List[Task]) -> List[Any]:
        """
        複数のタスクを並列実行します。

        Args:
            tasks: 実行するタスクのリスト

        Returns:
            タスクの実行結果のリスト
        """
        if not tasks:
            return []

        # Daskクライアントをセットアップ
        client = self._setup_client()

        try:
            # タスクを遅延実行オブジェクトに変換
            from dask.delayed import delayed

            print(f"並列計算を開始します（タスク数: {len(tasks)}）...")
            print(f"ワーカー数: {client.ncores}")

            # 各タスクをdelayedオブジェクトに変換
            delayed_results = []
            for i, task in enumerate(tasks):
                # タスクの状態を更新
                task.status = TaskStatus.RUNNING

                # 遅延実行オブジェクトを作成
                delayed_result = delayed(task.execute)()
                delayed_results.append(delayed_result)

            # 並列実行（非同期）
            print(f"タスクを並列実行中...")
            futures = client.compute(delayed_results)

            # 結果を取得（同期）- tqdmを使用して進捗状況と推定残り時間を表示
            print(f"結果を取得中...")
            # tqdmを使用して進捗バーを表示
            with tqdm(total=len(futures), desc="処理進捗", unit="タスク") as pbar:
                results = []
                for future in futures:
                    # 1つのタスクの結果を取得
                    result = client.gather(future)
                    results.append(result)
                    # 進捗バーを更新
                    pbar.update(1)

            # タスクの状態を更新
            for i, (task, result) in enumerate(zip(tasks, results)):
                task.status = TaskStatus.COMPLETED

            return results
        except Exception as e:
            # エラーハンドリング
            print(f"並列実行中にエラーが発生しました: {e}")
            # タスクの状態を更新
            for task in tasks:
                if task.status == TaskStatus.RUNNING:
                    task.status = TaskStatus.FAILED
            raise
        finally:
            # クライアントを閉じる
            client.close()

    def execute_async(self, task: Task) -> str:
        """
        タスクを非同期で実行します。

        Args:
            task: 実行するタスク

        Returns:
            非同期実行のジョブID
        """
        # 実装は後で行います
        return str(uuid.uuid4())  # ダミーの戻り値

    def get_execution_status(self, job_id: str) -> Dict[str, Any]:
        """
        非同期実行の状態を取得します。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            実行状態を表す辞書
        """
        # 実装は後で行います
        return {"status": "unknown", "job_id": job_id}  # ダミーの戻り値

    def get_execution_result(self, job_id: str) -> Any:
        """
        非同期実行の結果を取得します。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            タスクの実行結果
        """
        # 実装は後で行います
        return None  # ダミーの戻り値

    def cancel_execution(self, job_id: str) -> bool:
        """
        非同期実行をキャンセルします。

        Args:
            job_id: 非同期実行のジョブID

        Returns:
            キャンセルに成功した場合はTrue、そうでなければFalse
        """
        # 実装は後で行います
        return False  # ダミーの戻り値

    def _setup_local_cluster(self):
        """
        ローカル環境用のクラスタをセットアップします。

        Returns:
            設定されたローカルクラスタ
        """
        # 実装は後で行います
        pass

    def _setup_slurm_cluster(self):
        """
        Slurm環境用のクラスタをセットアップします。

        Returns:
            設定されたSlurmクラスタ
        """
        # 実装は後で行います
        pass

    def _setup_pbs_cluster(self):
        """
        PBS環境用のクラスタをセットアップします。

        Returns:
            設定されたPBSクラスタ
        """
        # 実装は後で行います
        pass

    def _optimize_performance(self, client):
        """
        Daskのパフォーマンスを最適化します。

        Args:
            client: Daskクライアント

        Returns:
            最適化されたDaskクライアント
        """
        # 実装は後で行います
        return client  # ダミーの戻り値
