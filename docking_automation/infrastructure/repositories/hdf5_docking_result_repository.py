import json  # metadataのシリアライズ/デシリアライズ用
import logging
import time  # リトライのためのsleep用
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

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
    上書きモードと追記モードをサポートしています。

    Attributes:
        hdf5_file_path (Path): HDF5ファイルのパス。
        lock_file_path (Path): ロックファイルのパス。
        mode (str): 保存モード。"overwrite"（上書き）または"append"（追記）。
    """

    def __init__(self, hdf5_file_path: Union[str, Path], mode: str = "overwrite") -> None:
        """リポジトリを初期化します。

        Args:
            hdf5_file_path (str | Path): HDF5ファイルのパス。
            mode (str, optional): 保存モード。"overwrite"（上書き）または"append"（追記）。
                デフォルトは"overwrite"。
        """
        self.hdf5_file_path = Path(hdf5_file_path)
        self.lock_file_path = self.hdf5_file_path.with_suffix(self.hdf5_file_path.suffix + ".lock")

        # モードの検証と設定
        if mode not in ["overwrite", "append"]:
            raise ValueError('モードは"overwrite"または"append"のいずれかを指定してください')
        self.mode = mode

        self._ensure_directory_exists()
        logger.info(
            f"HDF5リポジトリを初期化しました。ファイル: {self.hdf5_file_path}, "
            f"ロックファイル: {self.lock_file_path}, モード: {self.mode}"
        )

    def _ensure_directory_exists(self) -> None:
        """HDF5ファイルとロックファイルが存在するディレクトリを確認し、なければ作成します。"""
        self.hdf5_file_path.parent.mkdir(parents=True, exist_ok=True)
        self.lock_file_path.parent.mkdir(parents=True, exist_ok=True)  # 通常は同じディレクトリ

    def _exists(
        self,
        protein_content_hash: str,
        compound_content_hash: str,
    ) -> bool:
        """指定されたキーのデータが既に存在するかどうかを確認します。

        Args:
            protein_content_hash (str): タンパク質ファイルの内容ハッシュ値。
            compound_content_hash (str): 化合物のハッシュ値。

        Returns:
            bool: データが存在する場合はTrue、存在しない場合はFalse。

        Notes:
            TODO: protein_id, compound_set_id は現在使われていないので、関数の引数から削除可能
        """
        if not self.hdf5_file_path.exists():
            return False

        try:
            with h5py.File(self.hdf5_file_path, "r") as f:
                # 新しいパス形式
                new_group_path = f"/results/{protein_content_hash}/{compound_content_hash}"
                return new_group_path in f
        except Exception as e:
            logger.error(f"データの存在確認中にエラーが発生しました: {e}", exc_info=True)
            return False

    def save(self, docking_result: DockingResult) -> None:
        """ドッキング結果をHDF5ファイルに保存します。

        モードに応じて以下の動作をします：
        - "overwrite"モード: 既存の結果がある場合は上書きします。
        - "append"モード: 既存の結果がある場合はスキップします。

        ファイルロックを使用して排他制御を行い、ロック取得に失敗した場合はリトライします。

        Args:
            docking_result (DockingResult): 保存するドッキング結果。
        """
        lock = FileLock(self.lock_file_path)
        max_retries = 5  # 最大リトライ回数
        retry_delay = 1  # リトライ間隔（秒）

        # 追記モードで、既にデータが存在する場合はスキップ
        if self.mode == "append" and self._exists(
            docking_result.protein_content_hash,
            docking_result.compound_content_hash,
        ):
            logger.info(
                f"追記モード: 既存のデータが存在するためスキップします: "
                f"{docking_result.protein_content_hash}/{docking_result.compound_content_hash}"
            )
            return

        for attempt in range(max_retries):
            try:
                # タイムアウトを短めに設定し、リトライで対応
                with lock.acquire(timeout=10):
                    logger.debug(f"ロック取得成功 (試行 {attempt + 1}/{max_retries}): {self.lock_file_path}")
                    with h5py.File(self.hdf5_file_path, "a") as f:
                        # 新しいグループパス: /results/{protein_content_hash}/{compound_content_hash}
                        group_path = (
                            f"/results/{docking_result.protein_content_hash}/{docking_result.compound_content_hash}"
                        )

                        # 上書きモードの場合、既存のグループを削除
                        if self.mode == "overwrite" and group_path in f:
                            del f[group_path]

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
                        # SDFファイルの内容を読み込む
                        try:
                            with open(docking_result.result_path, "r") as f:
                                sdf_content = f.read()
                        except Exception as e:
                            logger.error(f"SDFファイルの読み込み中にエラーが発生しました: {e}", exc_info=True)
                            raise ValueError(f"SDFファイルの読み込みに失敗しました: {docking_result.result_path}")

                        # SDFファイルの内容を保存
                        if "sdf_content" in group:
                            del group["sdf_content"]
                        dt_str = h5py.string_dtype(encoding="utf-8")
                        group.create_dataset("sdf_content", data=sdf_content, dtype=dt_str)
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

                        logger.info(f"{self.mode}モード: ドッキング結果を保存しました: {group_path}")
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

    def load(self, result_id: str) -> Optional[DockingResult]:
        """
        指定されたIDに一致するドッキング結果をロードします。

        このメソッドは親クラスのインターフェースに準拠するために実装されています。
        result_idを解析して、protein_id、compound_set_id、compound_indexを抽出し、
        それらを使ってHDF5ファイルからデータを読み込みます。

        Args:
            result_id (str): 結果のID。形式は "{protein_id}_{compound_set_id}_{compound_index}" を想定。

        Returns:
            Optional[DockingResult]: 見つかったドッキング結果。見つからない場合はNone。
        """
        # result_idを解析
        parts = result_id.split("_")
        if len(parts) < 3:
            logger.warning(f"無効なresult_id形式です: {result_id}")
            return None

        # 最後の部分をcompound_indexとして解析
        try:
            compound_index = int(parts[-1])
        except ValueError:
            logger.warning(f"compound_indexを整数に変換できません: {parts[-1]}")
            return None

        # 残りの部分をprotein_idとcompound_set_idとして結合
        protein_id = parts[0]
        compound_set_id = "_".join(parts[1:-1])

        # HDF5ファイルからデータを読み込む
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return None

        lock = FileLock(self.lock_file_path)
        try:
            with lock.acquire(timeout=60):
                logger.debug(f"ロック取得 (読み込み): {self.lock_file_path}")
                with h5py.File(self.hdf5_file_path, "r") as f:
                    # 新しいパス形式で検索
                    found_group = None

                    if "results" in f:
                        results_group = f["results"]
                        # 全てのprotein_*グループを検索
                        for protein_group_name in results_group:
                            protein_group = results_group[protein_group_name]
                            # 全てのcompoundset_*グループを検索
                            for compound_set_group_name in protein_group:
                                compound_set_group = protein_group[compound_set_group_name]
                                # 指定されたインデックスを持つグループを検索
                                if str(compound_index) in compound_set_group:
                                    compound_group = compound_set_group[str(compound_index)]
                                    # 属性を確認
                                    if (
                                        compound_group.attrs["protein_id"] == protein_id
                                        and compound_group.attrs["compound_set_id"] == compound_set_id
                                    ):
                                        found_group = compound_group
                                        break

                            if found_group is not None:
                                break

                    # 新しいパス形式で検索（protein_content_hash/compound_content_hash）
                    if found_group is None:
                        # すべてのprotein_hash_*グループを検索
                        for protein_hash in f["results"]:
                            protein_group = f["results"][protein_hash]
                            # すべてのcompound_hash_*グループを検索
                            for compound_hash in protein_group:
                                compound_group = protein_group[compound_hash]
                                # 属性を確認
                                if (
                                    "protein_id" in compound_group.attrs
                                    and "compound_set_id" in compound_group.attrs
                                    and "compound_index" in compound_group.attrs
                                    and compound_group.attrs["protein_id"] == protein_id
                                    and compound_group.attrs["compound_set_id"] == compound_set_id
                                    and compound_group.attrs["compound_index"] == compound_index
                                ):
                                    found_group = compound_group
                                    break
                            if found_group is not None:
                                break

                        if found_group is None:
                            logger.debug(
                                f"指定された結果が見つかりません: {protein_id}_{compound_set_id}_{compound_index}"
                            )
                            return None

                    group = found_group

                    # データセットと属性から値を取得
                    docking_score = group["docking_score"][()]
                    sdf_content = group["sdf_content"][()].decode("utf-8")
                    metadata_json = group["metadata"][()].decode("utf-8")
                    metadata = json.loads(metadata_json)
                    result_id = group.attrs["id"]
                    version = group.attrs["version"]

                    # 一時ファイルを作成
                    import tempfile

                    with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
                        temp_file.write(sdf_content.encode("utf-8"))

                    # DockingResultオブジェクトを再構築
                    # protein_content_hashとcompoundset_content_hashを取得
                    # グループパスから抽出するか、デフォルト値を使用
                    group_path = group.name
                    parts = group_path.split("/")

                    # パスから抽出を試みる
                    protein_content_hash = "default_protein_hash"
                    compoundset_content_hash = "default_compound_hash"

                    for part in parts:
                        if part.startswith("protein_"):
                            protein_content_hash = part.replace("protein_", "")
                        elif part.startswith("compoundset_"):
                            compoundset_content_hash = part.replace("compoundset_", "")

                    # 化合物のハッシュ値を取得
                    compound_content_hash = "default_compound_hash"

                    # パスから抽出を試みる
                    parts = group_path.split("/")
                    if len(parts) >= 4:  # 新しいパス形式: /results/{protein_hash}/{compound_hash}
                        protein_content_hash = parts[2]
                        compound_content_hash = parts[3]

                    # 一時ファイルのパスを取得
                    temp_path = Path(temp_file.name)

                    result = DockingResult(
                        result_path=temp_path,
                        protein_id=group.attrs["protein_id"],
                        compound_set_id=group.attrs["compound_set_id"],
                        compound_index=group.attrs["compound_index"],
                        docking_score=docking_score,
                        protein_content_hash=protein_content_hash,
                        compound_content_hash=compound_content_hash,
                        compoundset_content_hash=compoundset_content_hash,
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
                    for protein_hash in results_group:
                        protein_group = results_group[protein_hash]
                        for compound_hash in protein_group:
                            compound_group = protein_group[compound_hash]

                            # データセットと属性から値を取得
                            # DockingResultの属性名に合わせて修正
                            docking_score = compound_group["docking_score"][()]
                            sdf_content = compound_group["sdf_content"][()].decode("utf-8")
                            metadata_json = compound_group["metadata"][()].decode("utf-8")
                            metadata = json.loads(metadata_json)
                            result_id = compound_group.attrs["id"]
                            version = compound_group.attrs["version"]

                            # 一時ファイルを作成
                            import tempfile

                            with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
                                temp_file.write(sdf_content.encode("utf-8"))

                            # DockingResultオブジェクトを再構築
                            # protein_content_hashとcompoundset_content_hashを取得
                            # グループパスから抽出するか、デフォルト値を使用
                            group_path = compound_group.name
                            parts = group_path.split("/")

                            # パスから抽出を試みる
                            protein_content_hash = parts[2]  # /results/{protein_hash}/{compound_hash}
                            compound_content_hash = parts[3]
                            compoundset_content_hash = "default_compound_hash"  # 後方互換性のため

                            # 一時ファイルのパスを取得
                            temp_path = Path(temp_file.name)

                            result = DockingResult(
                                result_path=temp_path,
                                protein_id=compound_group.attrs["protein_id"],
                                compound_set_id=compound_group.attrs["compound_set_id"],
                                compound_index=compound_group.attrs["compound_index"],  # 属性から取得
                                docking_score=docking_score,
                                protein_content_hash=protein_content_hash,
                                compound_content_hash=compound_content_hash,
                                compoundset_content_hash=compoundset_content_hash,
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
                    # 新しいパス形式で検索
                    found = False
                    group_path = ""

                    if "results" in f:
                        results_group = f["results"]
                        for protein_hash in results_group:
                            protein_group = results_group[protein_hash]
                            for compound_hash in protein_group:
                                compound_group = protein_group[compound_hash]
                                # 属性を確認
                                if (
                                    "protein_id" in compound_group.attrs
                                    and "compound_set_id" in compound_group.attrs
                                    and "compound_index" in compound_group.attrs
                                    and compound_group.attrs["protein_id"] == protein_id
                                    and compound_group.attrs["compound_set_id"] == compound_set_id
                                    and compound_group.attrs["compound_index"] == compound_index
                                ):
                                    group_path = compound_group.name
                                    found = True
                                    break
                            if found:
                                break

                    # 見つからなかった場合は従来のパス形式で検索（後方互換性のため）
                    if not found:
                        old_group_path = f"/results/{protein_id}/{compound_set_id}/{compound_index}"
                        if old_group_path in f:
                            group_path = old_group_path
                            found = True
                    if group_path in f:
                        del f[group_path]
                        logger.info(f"ドッキング結果を削除しました: {group_path}")
                    else:
                        logger.warning(f"削除対象の結果が見つかりません: {group_path}")

                    # 親グループが空になったら削除する (任意)
                    # 親グループのパスを取得
                    parent_path = "/".join(group_path.split("/")[:-1])
                    if parent_path in f and not list(f[parent_path].keys()):
                        del f[parent_path]
                        logger.debug(f"空の親グループを削除しました: {parent_path}")

                    # 祖父グループのパスを取得
                    grandparent_path = "/".join(parent_path.split("/")[:-1])
                    if grandparent_path in f and not list(f[grandparent_path].keys()):
                        del f[grandparent_path]
                        logger.debug(f"空の祖父グループを削除しました: {grandparent_path}")
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

        このメソッドは常に上書きモードで動作します（モード設定に関わらず）。

        Args:
            docking_result (DockingResult): 更新するドッキング結果。
        """
        # 現在のモードを一時的に保存
        current_mode = self.mode
        try:
            # 強制的に上書きモードに設定
            self.mode = "overwrite"
            logger.debug(
                f"updateメソッドが呼び出されました。強制的に上書きモードでsaveメソッドを実行します: "
                f"{docking_result.protein_id}/{docking_result.compound_set_id}/{docking_result.compound_index}"
            )
            self.save(docking_result)
        finally:
            # 元のモードに戻す
            self.mode = current_mode

    def get_repository_type(self) -> str:
        """リポジトリのタイプを返します。"""
        return "hdf5"

    def get_connection_details(self) -> Dict[str, Any]:
        """リポジトリへの接続詳細を返します。"""
        return {"file_path": str(self.hdf5_file_path)}
