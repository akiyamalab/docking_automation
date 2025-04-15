"""
Daskを使った並列処理Executorモジュール。

このモジュールは、Daskを使って複数のタスクを並列実行するためのExecutorを提供します。
ローカル環境とクラスタ環境（Slurm、PBS等）の両方に対応しています。
"""

import time
import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from dask.distributed import Client, LocalCluster, as_completed
from dask_jobqueue.pbs import PBSCluster
from dask_jobqueue.slurm import SLURMCluster
from tqdm import tqdm

from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
    DockingResultRepositoryFactory,
)

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

            cluster = LocalCluster(n_workers=self.n_workers, threads_per_worker=1)
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

    def execute_many(self, tasks: List[Task], **kwargs) -> List[Any]:
        """
        複数のタスクを並列実行します。

        Args:
            tasks: 実行するタスクのリスト
            **kwargs: 追加のパラメータ（例：repository_config）

        Returns:
            タスクの実行結果のリスト
        """
        # リポジトリ設定の検証
        repository_config = kwargs.get("repository_config", {})
        if not repository_config:
            raise ValueError("リポジトリ設定が指定されていません")
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

            # スケジューラ側で一元的に保存するためのリポジトリを作成
            repository_type = repository_config.get("repository_type")
            base_directory = repository_config.get("base_directory")
            config = repository_config.get("config", {})

            if not repository_type or not base_directory:
                raise ValueError("リポジトリ設定が不完全です。repository_typeとbase_directoryは必須です。")

            print(f"スケジューラ側で一元的に保存するためのリポジトリを作成します: {repository_type}")
            repository = DockingResultRepositoryFactory.create(
                repository_type=repository_type,
                base_directory=Path(base_directory),
                config=config,
            )

            # 結果を取得（同期）- tqdmを使用して進捗状況と推定残り時間を表示
            print(f"結果を取得中...")
            # tqdmを使用して進捗バーを表示
            with tqdm(total=len(futures), desc="処理進捗", unit="タスク") as pbar:
                results = [None] * len(futures)  # 結果を格納するリスト（順序を保持）
                future_to_idx = {f: i for i, f in enumerate(futures)}  # futureとインデックスのマッピング

                # 完了したタスクから順に結果を取得
                for future in as_completed(futures):
                    idx = future_to_idx[future]
                    task_results = client.gather(future)
                    results[idx] = task_results

                    # スケジューラ側で一元的に保存
                    # DockingResultCollectionの場合は各結果を個別に保存
                    if hasattr(task_results, "__iter__") and not isinstance(task_results, str):
                        for result in task_results:
                            try:
                                repository.save(result)
                            except Exception as e:
                                print(
                                    f"警告: 結果の保存中にエラーが発生しました - "
                                    f"タンパク質ID={result.protein_id if hasattr(result, 'protein_id') else 'unknown'}, "
                                    f"化合物セットID={result.compound_set_id if hasattr(result, 'compound_set_id') else 'unknown'}, "
                                    f"化合物インデックス={result.compound_index if hasattr(result, 'compound_index') else 'unknown'}, "
                                    f"エラー: {e}"
                                )

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
            # クライアントを閉じる前に少し待機（ワーカーのハートビート通信が完了するのを待つ）
            time.sleep(1)
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
