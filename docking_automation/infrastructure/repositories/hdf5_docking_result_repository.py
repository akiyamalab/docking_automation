import json  # metadataのシリアライズ/デシリアライズ用
import logging
import time  # リトライのためのsleep用
from pathlib import Path
from typing import Any, Dict, List, Optional, cast

import h5py
from filelock import FileLock, Timeout  # Timeout例外をインポート

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepository,
)

logger = logging.getLogger(__name__)


class HDF5DockingResultRepository(DockingResultRepository):
    """HDF5ファイルを使用してドッキング結果を永続化するリポジトリ。

    ファイルロックを使用して、複数のプロセスからの同時書き込みを制御します。

    Attributes:
        hdf5_file_path (Path): HDF5ファイルのパス。
        lock_file_path (Path): ロックファイルのパス。
    """

    def __init__(self, hdf5_file_path: str | Path):
        """リポジトリを初期化します。

        Args:
            hdf5_file_path (str | Path): HDF5ファイルのパス。
        """
        self.hdf5_file_path = Path(hdf5_file_path)
        self.lock_file_path = self.hdf5_file_path.with_suffix(self.hdf5_file_path.suffix + ".lock")
        self._ensure_directory_exists()
        logger.info(
            f"HDF5リポジトリを初期化しました。ファイル: {self.hdf5_file_path}, ロックファイル: {self.lock_file_path}"
        )

    def _ensure_directory_exists(self):
        """HDF5ファイルとロックファイルが存在するディレクトリを確認し、なければ作成します。"""
        self.hdf5_file_path.parent.mkdir(parents=True, exist_ok=True)
        self.lock_file_path.parent.mkdir(parents=True, exist_ok=True)  # 通常は同じディレクトリ

    def save(self, docking_result: DockingResult) -> None:
        """単一のドッキング結果をHDF5ファイルに保存します。

        既存の結果がある場合は上書きされます。ファイルロックを使用して排他制御を行い、
        ロック取得に失敗した場合はリトライします。

        Args:
            docking_result (DockingResult): 保存するドッキング結果。
        """
        lock = FileLock(self.lock_file_path)
        max_retries = 5  # 最大リトライ回数
        retry_delay = 1  # リトライ間隔（秒）

        for attempt in range(max_retries):
            try:
                # タイムアウトを短めに設定し、リトライで対応
                with lock.acquire(timeout=10):
                    logger.debug(f"ロック取得成功 (試行 {attempt + 1}/{max_retries}): {self.lock_file_path}")
                    with h5py.File(self.hdf5_file_path, "a") as f:
                        # グループパス: /results/{protein_id}/{compound_set_id}/{compound_index}
                        group_path = f"/results/{docking_result.protein_id}/{docking_result.compound_set_id}/{docking_result.compound_index}"
                        group = f.require_group(group_path)

                        # 属性やデータセットとして情報を保存 (既存データの上書き)
                        if "protein_id" in group.attrs:
                            del group.attrs["protein_id"]
                        group.attrs["protein_id"] = docking_result.protein_id
                        if "compound_set_id" in group.attrs:
                            del group.attrs["compound_set_id"]
                        group.attrs["compound_set_id"] = docking_result.compound_set_id
                        if "compound_index" in group.attrs:
                            del group.attrs["compound_index"]
                        group.attrs["compound_index"] = docking_result.compound_index
                        if "docking_score" in group:
                            del group["docking_score"]
                        group.create_dataset("docking_score", data=docking_result.docking_score, dtype="f8")
                        if "result_path" in group:
                            del group["result_path"]
                        dt_str = h5py.string_dtype(encoding="utf-8")
                        group.create_dataset("result_path", data=str(docking_result.result_path), dtype=dt_str)
                        if "metadata" in group:
                            del group["metadata"]
                        metadata_json = json.dumps(docking_result.metadata)
                        group.create_dataset("metadata", data=metadata_json, dtype=dt_str)
                        if "id" in group.attrs:
                            del group.attrs["id"]
                        group.attrs["id"] = docking_result.id
                        if "version" in group.attrs:
                            del group.attrs["version"]
                        group.attrs["version"] = docking_result.version

                        logger.info(f"ドッキング結果を保存しました: {group_path}")
                        return  # 保存成功したらループを抜ける

            except Timeout:  # filelock.Timeout 例外をキャッチ
                logger.warning(
                    f"ロック取得タイムアウト (試行 {attempt + 1}/{max_retries})。{retry_delay}秒待機してリトライします: {self.lock_file_path}"
                )
                if attempt < max_retries - 1:
                    time.sleep(retry_delay + (attempt * 0.5))  # リトライごとに少し待機時間を増やす
                else:
                    logger.error(f"ロック取得に失敗しました ({max_retries}回リトライ後): {self.lock_file_path}")
                    raise  # 最大リトライ回数を超えたらエラーを送出
            except Exception as e:
                logger.error(f"HDF5ファイルへの保存中に予期せぬエラーが発生しました: {e}", exc_info=True)
                raise  # その他のエラーはそのまま送出
            finally:
                # finallyブロックはtry-exceptブロック全体に対して実行されるため、
                # ロックが取得されている場合のみ解放する
                if lock.is_locked:
                    lock.release()
                    logger.debug(f"ロック解放 (試行 {attempt + 1}): {self.lock_file_path}")

    def load(self, protein_id: str, compound_set_id: str, compound_index: int) -> Optional[DockingResult]:
        """指定されたIDとインデックスに一致するドッキング結果をロードします。

        ファイルロックを使用して読み込み中の書き込みを防ぎます。

        Args:
            protein_id (str): タンパク質ID。
            compound_set_id (str): 化合物セットID。
            compound_index (int): 化合物インデックス。

        Returns:
            Optional[DockingResult]: 見つかったドッキング結果。見つからない場合はNone。
        """
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return None

        lock = FileLock(self.lock_file_path)
        try:
            with lock.acquire(timeout=60):
                logger.debug(f"ロック取得 (読み込み): {self.lock_file_path}")
                with h5py.File(self.hdf5_file_path, "r") as f:
                    # グループパスを修正
                    group_path = f"/results/{protein_id}/{compound_set_id}/{compound_index}"
                    if group_path not in f:
                        logger.debug(f"指定された結果が見つかりません: {group_path}")
                        return None

                    group = f[group_path]

                    # データセットと属性から値を取得
                    # DockingResultの属性名に合わせて修正
                    docking_score = group["docking_score"][()]
                    result_path_str = group["result_path"][()].decode("utf-8")
                    metadata_json = group["metadata"][()].decode("utf-8")
                    metadata = json.loads(metadata_json)
                    result_id = group.attrs["id"]
                    version = group.attrs["version"]

                    # DockingResultオブジェクトを再構築 (コンストラクタを使用)
                    # result_pathをPathオブジェクトに戻す
                    result = DockingResult(
                        result_path=Path(result_path_str),
                        protein_id=group.attrs["protein_id"],
                        compound_set_id=group.attrs["compound_set_id"],
                        compound_index=group.attrs["compound_index"],
                        docking_score=docking_score,
                        metadata=metadata,
                        id=result_id,
                        version=version,
                    )
                    logger.info(f"ドッキング結果をロードしました: {group_path}")
                    return result

        except TimeoutError:
            logger.error(f"ロックの取得にタイムアウトしました (読み込み): {self.lock_file_path}")
            raise
        except Exception as e:
            logger.error(f"HDF5ファイルからの読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise
        finally:
            if lock.is_locked:
                lock.release()
                logger.debug(f"ロック解放 (読み込み): {self.lock_file_path}")

    def load_all(self) -> DockingResultCollection:
        """HDF5ファイルに保存されているすべてのドッキング結果をロードします。

        Returns:
            DockingResultCollection: すべてのドッキング結果を含むコレクション。
        """
        results: List[DockingResult] = []
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return DockingResultCollection(results)

        lock = FileLock(self.lock_file_path)
        try:
            with lock.acquire(timeout=60):
                logger.debug(f"ロック取得 (全件読み込み): {self.lock_file_path}")
                with h5py.File(self.hdf5_file_path, "r") as f:
                    if "results" not in f:
                        logger.info("結果グループが存在しません。空のコレクションを返します。")
                        return DockingResultCollection(results)

                    results_group = f["results"]
                    for protein_id in results_group:
                        protein_group = results_group[protein_id]
                        for compound_set_id in protein_group:
                            compound_set_group = protein_group[compound_set_id]
                            for compound_index_str in compound_set_group:  # インデックスは文字列キーになっている
                                compound_group = compound_set_group[compound_index_str]

                                # データセットと属性から値を取得
                                # DockingResultの属性名に合わせて修正
                                docking_score = compound_group["docking_score"][()]
                                result_path_str = compound_group["result_path"][()].decode("utf-8")
                                metadata_json = compound_group["metadata"][()].decode("utf-8")
                                metadata = json.loads(metadata_json)
                                result_id = compound_group.attrs["id"]
                                version = compound_group.attrs["version"]

                                # DockingResultオブジェクトを再構築
                                result = DockingResult(
                                    result_path=Path(result_path_str),
                                    protein_id=compound_group.attrs["protein_id"],
                                    compound_set_id=compound_group.attrs["compound_set_id"],
                                    compound_index=compound_group.attrs["compound_index"],  # 属性から取得
                                    docking_score=docking_score,
                                    metadata=metadata,
                                    id=result_id,
                                    version=version,
                                )
                                results.append(result)
                    logger.info(f"{len(results)}件のドッキング結果をロードしました。")

        except TimeoutError:
            logger.error(f"ロックの取得にタイムアウトしました (全件読み込み): {self.lock_file_path}")
            raise
        except Exception as e:
            logger.error(f"HDF5ファイルからの全件読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise
        finally:
            if lock.is_locked:
                lock.release()
                logger.debug(f"ロック解放 (全件読み込み): {self.lock_file_path}")

        return DockingResultCollection(results)

    def delete(self, protein_id: str, compound_set_id: str, compound_index: int) -> None:
        """指定されたドッキング結果を削除します。

        Args:
            protein_id (str): タンパク質ID。
            compound_set_id (str): 化合物セットID。
            compound_index (int): 化合物インデックス。
        """
        lock = FileLock(self.lock_file_path)
        try:
            with lock.acquire(timeout=60):
                logger.debug(f"ロック取得 (削除): {self.lock_file_path}")
                with h5py.File(self.hdf5_file_path, "a") as f:  # 書き込みモードで開く
                    # グループパスを修正
                    group_path = f"/results/{protein_id}/{compound_set_id}/{compound_index}"
                    if group_path in f:
                        del f[group_path]
                        logger.info(f"ドッキング結果を削除しました: {group_path}")
                    else:
                        logger.warning(f"削除対象の結果が見つかりません: {group_path}")

                    # 親グループが空になったら削除する (任意)
                    compound_set_group_path = f"/results/{protein_id}/{compound_set_id}"
                    if compound_set_group_path in f and not list(f[compound_set_group_path].keys()):
                        del f[compound_set_group_path]
                        logger.debug(f"空の化合物セットグループを削除しました: {compound_set_group_path}")

                    protein_group_path = f"/results/{protein_id}"
                    if protein_group_path in f and not list(f[protein_group_path].keys()):
                        del f[protein_group_path]
                        logger.debug(f"空のタパク質グループを削除しました: {protein_group_path}")
                    if "/results" in f and not list(f["/results"].keys()):
                        del f["/results"]
                        logger.debug("空の結果グループを削除しました: /results")

        except TimeoutError:
            logger.error(f"ロックの取得にタイムアウトしました (削除): {self.lock_file_path}")
            raise
        except Exception as e:
            logger.error(f"HDF5ファイルからの削除中にエラーが発生しました: {e}", exc_info=True)
            raise
        finally:
            if lock.is_locked:
                lock.release()
                logger.debug(f"ロック解放 (削除): {self.lock_file_path}")

    def update(self, docking_result: DockingResult) -> None:
        """既存のドッキング結果を更新します。

        saveメソッドと同じ動作（上書き）になります。

        Args:
            docking_result (DockingResult): 更新するドッキング結果。
        """
        # ログ出力のキーを修正
        logger.debug(
            f"updateメソッドが呼び出されました。saveメソッドを実行します: {docking_result.protein_id}/{docking_result.compound_set_id}/{docking_result.compound_index}"
        )
        self.save(docking_result)

    def get_repository_type(self) -> str:
        """リポジトリのタイプを返します。"""
        return "hdf5"

    def get_connection_details(self) -> Dict[str, Any]:
        """リポジトリへの接続詳細を返します。"""
        return {"file_path": str(self.hdf5_file_path)}
