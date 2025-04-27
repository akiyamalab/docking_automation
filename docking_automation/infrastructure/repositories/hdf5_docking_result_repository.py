import json  # metadataのシリアライズ/デシリアライズ用
import logging
from pathlib import Path
from typing import Any, Dict, List, Optional, Union

import h5py

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepository,
)

logger = logging.getLogger(__name__)


class HDF5DockingResultRepository(DockingResultRepository):
    """HDF5ファイルを使用してドッキング結果を永続化するリポジトリ。

    SWMRモード（Single Writer Multiple Reader）を使用して、複数のプロセスからの同時アクセスを制御します。
    上書きモードと追記モードをサポートしています。

    Attributes:
        hdf5_file_path (Path): HDF5ファイルのパス。
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

        # モードの検証と設定
        if mode not in ["overwrite", "append"]:
            raise ValueError('モードは"overwrite"または"append"のいずれかを指定してください')
        self.mode = mode

        self._ensure_directory_exists()
        logger.info(f"HDF5リポジトリを初期化しました。ファイル: {self.hdf5_file_path}, モード: {self.mode}")

    def _ensure_directory_exists(self) -> None:
        """HDF5ファイルが存在するディレクトリを確認し、なければ作成します。"""
        self.hdf5_file_path.parent.mkdir(parents=True, exist_ok=True)

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
            # swmrモード（Single Writer Multiple Reader）を使用
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
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

        SWMRモードを使用して、複数のプロセスからの同時アクセスを制御します。

        Args:
            docking_result (DockingResult): 保存するドッキング結果。
        """
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

        try:
            # SWMRモードで書き込み
            with h5py.File(self.hdf5_file_path, "a", libver="latest") as f:
                f.swmr_mode = True
                # 新しいグループパス: /results/{protein_content_hash}/{compound_content_hash}
                group_path = f"/results/{docking_result.protein_content_hash}/{docking_result.compound_content_hash}"

                # 上書きモードの場合、既存のグループを削除
                if self.mode == "overwrite" and group_path in f:
                    del f[group_path]

                group = f.require_group(group_path)

                # 属性やデータセットとして情報を保存
                group.attrs["protein_id"] = docking_result.protein_id
                group.attrs["compound_set_id"] = docking_result.compound_set_id
                group.attrs["compound_index"] = docking_result.compound_index
                group.create_dataset("docking_score", data=docking_result.docking_score, dtype="f8")

                # SDFファイルの内容を読み込む
                try:
                    with open(docking_result.result_path, "r") as sdf_file:
                        sdf_content = sdf_file.read()
                except Exception as e:
                    logger.error(f"SDFファイルの読み込み中にエラーが発生しました: {e}", exc_info=True)
                    raise ValueError(f"SDFファイルの読み込みに失敗しました: {docking_result.result_path}")

                # SDFファイルの内容を保存
                dt_str = h5py.string_dtype(encoding="utf-8")
                group.create_dataset("sdf_content", data=sdf_content, dtype=dt_str)
                metadata_json = json.dumps(docking_result.metadata)
                group.create_dataset("metadata", data=metadata_json, dtype=dt_str)
                group.attrs["id"] = docking_result.id
                group.attrs["version"] = docking_result.version

                logger.info(f"{self.mode}モード: ドッキング結果を保存しました: {group_path}")
        except Exception as e:
            logger.error(f"HDF5ファイルへの保存中に予期せぬエラーが発生しました: {e}", exc_info=True)
            raise

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

        try:
            # SWMRモードで読み込み
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
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
                        logger.debug(f"指定された結果が見つかりません: {protein_id}_{compound_set_id}_{compound_index}")
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

        except Exception as e:
            logger.error(f"HDF5ファイルからの読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

    def load_by_hashes(self, protein_content_hash: str, compound_content_hash: str) -> Optional[DockingResult]:
        """タンパク質と化合物のハッシュ値に基づいてドッキング結果をロードします。

        Args:
            protein_content_hash (str): タンパク質ファイルの内容ハッシュ値。
            compound_content_hash (str): 化合物のハッシュ値。

        Returns:
            Optional[DockingResult]: 見つかったドッキング結果。見つからない場合はNone。
        """
        # ファイルが存在するか確認
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return None

        # 指定されたハッシュ値のデータが存在するか確認
        if not self._exists(protein_content_hash, compound_content_hash):
            logger.debug(f"指定されたハッシュ値のデータが存在しません: {protein_content_hash}/{compound_content_hash}")
            return None

        try:
            # SWMRモードで読み込み
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                # 新しいパス形式で検索
                group_path = f"/results/{protein_content_hash}/{compound_content_hash}"

                if group_path not in f:
                    logger.debug(f"指定されたパスが存在しません: {group_path}")
                    return None

                group = f[group_path]

                # データセットと属性から値を取得
                docking_score = group["docking_score"][()]
                sdf_content = group["sdf_content"][()].decode("utf-8")
                metadata_json = group["metadata"][()].decode("utf-8")
                metadata = json.loads(metadata_json)
                result_id = group.attrs["id"]
                version = group.attrs["version"]
                protein_id = group.attrs["protein_id"]
                compound_set_id = group.attrs["compound_set_id"]
                compound_index = group.attrs["compound_index"]

                # 一時ファイルを作成
                import tempfile

                with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
                    temp_file.write(sdf_content.encode("utf-8"))

                # 一時ファイルのパスを取得
                temp_path = Path(temp_file.name)

                # 再利用されたことを示すメタデータを追加
                metadata["reused"] = True

                result = DockingResult(
                    result_path=temp_path,
                    protein_id=protein_id,
                    compound_set_id=compound_set_id,
                    compound_index=compound_index,
                    docking_score=docking_score,
                    protein_content_hash=protein_content_hash,
                    compound_content_hash=compound_content_hash,
                    compoundset_content_hash=metadata.get("compoundset_content_hash", ""),
                    metadata=metadata,
                    id=result_id,
                    version=version,
                )
                logger.info(f"ハッシュ値によるドッキング結果のロードに成功しました: {group_path}")
                return result

        except Exception as e:
            logger.error(f"HDF5ファイルからのハッシュ値による読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

    def load_all(self) -> DockingResultCollection:
        """HDF5ファイルに保存されているすべてのドッキング結果をロードします。

        Returns:
            DockingResultCollection: すべてのドッキング結果を含むコレクション。
        """
        results: List[DockingResult] = []
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return DockingResultCollection(results)

        try:
            # SWMRモードで読み込み
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
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

        except Exception as e:
            logger.error(f"HDF5ファイルからの全件読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

        return DockingResultCollection(results)

    def delete(self, protein_id: str, compound_set_id: str, compound_index: int) -> None:
        """指定されたドッキング結果を削除します。

        Args:
            protein_id (str): タンパク質ID。
            compound_set_id (str): 化合物セットID。
            compound_index (int): 化合物インデックス。
        """
        try:
            # SWMRモードで書き込み
            with h5py.File(self.hdf5_file_path, "a", libver="latest") as f:
                f.swmr_mode = True
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

        except Exception as e:
            logger.error(f"HDF5ファイルからの削除中にエラーが発生しました: {e}", exc_info=True)
            raise

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
