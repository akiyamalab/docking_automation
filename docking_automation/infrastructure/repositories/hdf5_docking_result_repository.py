import gzip
import logging
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set, Tuple, Union

import h5py
import numpy as np

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepository,
)

logger = logging.getLogger(__name__)

_FAILED_SCORE_SENTINEL = -999.0


class HDF5DockingResultRepository(DockingResultRepository):
    """HDF5ファイルを使用してドッキング結果を永続化するリポジトリ。

    SWMRモード（Single Writer Multiple Reader）を使用して、複数のプロセスからの同時アクセスを制御します。
    上書きモードと追記モードをサポートしています。

    Attributes:
        hdf5_file_path (Path): HDF5ファイルのパス。
        mode (str): 保存モード。"overwrite"（上書き）または"append"（追記）。
        schema_version (str): "v2"=per-pair スキーマ（後方互換）, "v3"=protein-bundle スキーマ。
    """

    def __init__(
        self,
        hdf5_file_path: Union[str, Path],
        mode: str = "overwrite",
        schema_version: str = "v2",
    ) -> None:
        self.hdf5_file_path = Path(hdf5_file_path)

        if mode not in ["overwrite", "append"]:
            raise ValueError('モードは"overwrite"または"append"のいずれかを指定してください')
        self.mode = mode

        if schema_version not in ["v2", "v3"]:
            raise ValueError('schema_versionは"v2"または"v3"のいずれかを指定してください')
        self.schema_version = schema_version

        self._ensure_directory_exists()
        logger.info(
            f"HDF5リポジトリを初期化しました。ファイル: {self.hdf5_file_path}, "
            f"モード: {self.mode}, スキーマ: {self.schema_version}"
        )

    def _ensure_directory_exists(self) -> None:
        self.hdf5_file_path.parent.mkdir(parents=True, exist_ok=True)

    def _exists(
        self,
        protein_content_hash: str,
        compound_content_hash: str,
    ) -> bool:
        if not self.hdf5_file_path.exists():
            return False

        try:
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                group_path = f"/results/{protein_content_hash}/{compound_content_hash}"
                return group_path in f
        except Exception as e:
            logger.error(f"データの存在確認中にエラーが発生しました: {e}", exc_info=True)
            return False

    def save(self, docking_result: DockingResult) -> None:
        """ドッキング結果をHDF5ファイルに保存します (Phase 2スキーマ)。

        pose_blob: gzip level=9 圧縮済みSDF bytes
        docking_score: float32
        computed_at: ISO8601 UTC文字列
        source: "vina" (固定)
        top_n: int32(1) (固定)
        """
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
            with h5py.File(self.hdf5_file_path, "a", libver="latest") as f:
                f.swmr_mode = True
                group_path = f"/results/{docking_result.protein_content_hash}/{docking_result.compound_content_hash}"

                if self.mode == "overwrite" and group_path in f:
                    del f[group_path]

                group = f.require_group(group_path)

                group.attrs["protein_id"] = docking_result.protein_id
                group.attrs["compound_set_id"] = docking_result.compound_set_id
                group.attrs["compound_index"] = docking_result.compound_index
                group.attrs["id"] = docking_result.id
                group.attrs["version"] = docking_result.version

                try:
                    with open(docking_result.result_path, "rb") as sdf_file:
                        sdf_bytes = sdf_file.read()
                except Exception as e:
                    logger.error(f"SDFファイルの読み込み中にエラーが発生しました: {e}", exc_info=True)
                    raise ValueError(f"SDFファイルの読み込みに失敗しました: {docking_result.result_path}")

                pose_blob = gzip.compress(sdf_bytes, compresslevel=9)

                group.create_dataset("pose_blob", data=np.frombuffer(pose_blob, dtype=np.uint8))
                group.create_dataset("docking_score", data=np.float32(docking_result.docking_score))
                dt_str = h5py.string_dtype(encoding="utf-8")
                group.create_dataset("computed_at", data=datetime.utcnow().isoformat() + "Z", dtype=dt_str)
                group.create_dataset("source", data="vina", dtype=dt_str)
                group.create_dataset("top_n", data=np.int32(1))

                logger.info(f"{self.mode}モード: ドッキング結果を保存しました: {group_path}")
        except Exception as e:
            logger.error(f"HDF5ファイルへの保存中に予期せぬエラーが発生しました: {e}", exc_info=True)
            raise

    def _decode_group(self, group: Any, protein_content_hash: str, compound_content_hash: str) -> DockingResult:
        """HDF5グループからDockingResultを復元する共通ロジック。

        Phase 2スキーマ（pose_blob）と旧スキーマ（sdf_content）の両方に対応。
        """
        docking_score = float(group["docking_score"][()])

        if "pose_blob" in group:
            pose_blob_bytes = group["pose_blob"][()].tobytes()
            sdf_content = gzip.decompress(pose_blob_bytes).decode("utf-8")
        elif "sdf_content" in group:
            sdf_content = group["sdf_content"][()].decode("utf-8")
        else:
            sdf_content = ""

        result_id = group.attrs.get("id", f"{protein_content_hash}_{compound_content_hash}")
        version = int(group.attrs.get("version", 1))
        protein_id = group.attrs.get("protein_id", "")
        compound_set_id = group.attrs.get("compound_set_id", "")
        compound_index = int(group.attrs.get("compound_index", 0))

        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False, mode="w", encoding="utf-8") as temp_file:
            temp_file.write(sdf_content)
            temp_path = Path(temp_file.name)

        return DockingResult(
            result_path=temp_path,
            protein_id=protein_id,
            compound_set_id=compound_set_id,
            compound_index=compound_index,
            docking_score=docking_score,
            protein_content_hash=protein_content_hash,
            compound_content_hash=compound_content_hash,
            compoundset_content_hash=compound_content_hash,
            metadata={},
            id=result_id,
            version=version,
        )

    def load(self, result_id: str) -> Optional[DockingResult]:
        """指定されたIDに一致するドッキング結果をロードします。

        result_idを解析して、protein_id、compound_set_id、compound_indexを抽出し、
        それらを使ってHDF5ファイルからデータを読み込みます。

        Args:
            result_id (str): 結果のID。形式は "{protein_id}_{compound_set_id}_{compound_index}" を想定。
        """
        parts = result_id.split("_")
        if len(parts) < 3:
            logger.warning(f"無効なresult_id形式です: {result_id}")
            return None

        try:
            compound_index = int(parts[-1])
        except ValueError:
            logger.warning(f"compound_indexを整数に変換できません: {parts[-1]}")
            return None

        protein_id = parts[0]
        compound_set_id = "_".join(parts[1:-1])

        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return None

        try:
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                if "results" not in f:
                    return None

                found_group = None
                found_protein_hash = ""
                found_compound_hash = ""

                for protein_hash in f["results"]:
                    protein_group = f["results"][protein_hash]
                    for compound_hash in protein_group:
                        compound_group = protein_group[compound_hash]
                        if (
                            "protein_id" in compound_group.attrs
                            and "compound_set_id" in compound_group.attrs
                            and "compound_index" in compound_group.attrs
                            and compound_group.attrs["protein_id"] == protein_id
                            and compound_group.attrs["compound_set_id"] == compound_set_id
                            and int(compound_group.attrs["compound_index"]) == compound_index
                        ):
                            found_group = compound_group
                            found_protein_hash = protein_hash
                            found_compound_hash = compound_hash
                            break
                    if found_group is not None:
                        break

                if found_group is None:
                    logger.debug(f"指定された結果が見つかりません: {result_id}")
                    return None

                result = self._decode_group(found_group, found_protein_hash, found_compound_hash)
                logger.info(f"ドッキング結果をロードしました: {found_group.name}")
                return result

        except Exception as e:
            logger.error(f"HDF5ファイルからの読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

    def load_by_hashes(self, protein_content_hash: str, compound_content_hash: str) -> Optional[DockingResult]:
        """タンパク質と化合物のハッシュ値に基づいてドッキング結果をロードします。"""
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return None

        if not self._exists(protein_content_hash, compound_content_hash):
            logger.debug(f"指定されたハッシュ値のデータが存在しません: {protein_content_hash}/{compound_content_hash}")
            return None

        try:
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                group_path = f"/results/{protein_content_hash}/{compound_content_hash}"
                if group_path not in f:
                    return None

                group = f[group_path]
                result = self._decode_group(group, protein_content_hash, compound_content_hash)
                logger.info(f"ハッシュ値によるドッキング結果のロードに成功しました: {group_path}")
                return result

        except Exception as e:
            logger.error(f"HDF5ファイルからのハッシュ値による読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

    def load_all(self) -> DockingResultCollection:
        """HDF5ファイルに保存されているすべてのドッキング結果をロードします。"""
        results: List[DockingResult] = []
        if not self.hdf5_file_path.exists():
            logger.warning(f"HDF5ファイルが存在しません: {self.hdf5_file_path}")
            return DockingResultCollection(results)

        try:
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                if "results" not in f:
                    logger.info("結果グループが存在しません。空のコレクションを返します。")
                    return DockingResultCollection(results)

                for protein_hash in f["results"]:
                    protein_group = f["results"][protein_hash]
                    for compound_hash in protein_group:
                        compound_group = protein_group[compound_hash]
                        result = self._decode_group(compound_group, protein_hash, compound_hash)
                        results.append(result)

                logger.info(f"{len(results)}件のドッキング結果をロードしました。")

        except Exception as e:
            logger.error(f"HDF5ファイルからの全件読み込み中にエラーが発生しました: {e}", exc_info=True)
            raise

        return DockingResultCollection(results)

    def get_all_keys(self) -> Set[Tuple[str, str]]:
        """全既存 (protein_hash, compound_hash) ペアを返す。ScreeningRunnerの差分検出に使用。"""
        keys: Set[Tuple[str, str]] = set()
        if not self.hdf5_file_path.exists():
            return keys

        try:
            with h5py.File(self.hdf5_file_path, "r", swmr=True) as f:
                if "results" not in f:
                    return keys
                for protein_hash in f["results"]:
                    for compound_hash in f["results"][protein_hash]:
                        keys.add((protein_hash, compound_hash))
        except Exception as e:
            logger.error(f"get_all_keys中にエラーが発生しました: {e}", exc_info=True)
            raise

        return keys

    def delete(self, protein_id: str, compound_set_id: str, compound_index: int) -> None:
        """指定されたドッキング結果を削除します。"""
        try:
            with h5py.File(self.hdf5_file_path, "a", libver="latest") as f:
                f.swmr_mode = True
                found = False
                group_path = ""

                if "results" in f:
                    results_group = f["results"]
                    for protein_hash in results_group:
                        protein_group = results_group[protein_hash]
                        for compound_hash in protein_group:
                            compound_group = protein_group[compound_hash]
                            if (
                                "protein_id" in compound_group.attrs
                                and "compound_set_id" in compound_group.attrs
                                and "compound_index" in compound_group.attrs
                                and compound_group.attrs["protein_id"] == protein_id
                                and compound_group.attrs["compound_set_id"] == compound_set_id
                                and int(compound_group.attrs["compound_index"]) == compound_index
                            ):
                                group_path = compound_group.name
                                found = True
                                break
                        if found:
                            break

                if group_path and group_path in f:
                    del f[group_path]
                    logger.info(f"ドッキング結果を削除しました: {group_path}")
                else:
                    logger.warning(f"削除対象の結果が見つかりません: {protein_id}/{compound_set_id}/{compound_index}")

                if group_path:
                    parent_path = "/".join(group_path.split("/")[:-1])
                    if parent_path in f and not list(f[parent_path].keys()):
                        del f[parent_path]
                        logger.debug(f"空の親グループを削除しました: {parent_path}")

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
        """既存のドッキング結果を更新します（常に上書きモードで動作）。"""
        current_mode = self.mode
        try:
            self.mode = "overwrite"
            self.save(docking_result)
        finally:
            self.mode = current_mode

    def get_repository_type(self) -> str:
        return "hdf5"

    def get_connection_details(self) -> Dict[str, Any]:
        return {"file_path": str(self.hdf5_file_path)}

    # ── v3 protein-bundle スキーマ ────────────────────────────────────────────

    def write_bundle(self, protein_hash: str, entries: List[Dict[str, Any]]) -> int:
        """protein-bundle スキーマ (v3) で 1 タンパク質分の結果を書き込む。

        entries: list of dicts with keys:
          - compound_hash: str
          - score: Optional[float]  (None → _FAILED_SCORE_SENTINEL = -999.0)
          - pose_blob: Optional[bytes]  (None → b"")
          - source: str  (default "vina")
          - top_n: int   (default 1)

        Returns: number of newly written entries (skips duplicates in append mode).
        """
        if not entries:
            return 0

        group_path = f"/results/{protein_hash}"
        dt_str = h5py.string_dtype(encoding="utf-8")
        vlen_bytes = h5py.vlen_dtype(np.uint8)

        with h5py.File(self.hdf5_file_path, "a", libver="latest") as f:
            if group_path in f:
                grp = f[group_path]
                if "compound_hashes" in grp:
                    existing_hashes: Set[str] = set(
                        h.decode("utf-8") if isinstance(h, bytes) else str(h)
                        for h in grp["compound_hashes"][:]
                    )
                else:
                    existing_hashes = set()
            else:
                grp = f.require_group(group_path)
                existing_hashes = set()

            new_entries = [e for e in entries if e["compound_hash"] not in existing_hashes]
            if not new_entries:
                return 0

            compound_hashes = [e["compound_hash"] for e in new_entries]
            scores = np.array(
                [
                    e["score"] if e.get("score") is not None else _FAILED_SCORE_SENTINEL
                    for e in new_entries
                ],
                dtype=np.float32,
            )
            pose_blobs_raw = [e.get("pose_blob") or b"" for e in new_entries]
            computed_ats = [datetime.utcnow().isoformat() + "Z" for _ in new_entries]
            sources = [e.get("source", "vina") for e in new_entries]
            top_ns = np.array([e.get("top_n", 1) for e in new_entries], dtype=np.int32)
            n_new = len(new_entries)

            if "compound_hashes" in grp:
                n_exist = grp["compound_hashes"].shape[0]
                n_total = n_exist + n_new

                grp["compound_hashes"].resize(n_total, axis=0)
                grp["compound_hashes"][n_exist:] = [c.encode("utf-8") for c in compound_hashes]

                grp["docking_scores"].resize(n_total, axis=0)
                grp["docking_scores"][n_exist:] = scores

                grp["pose_blobs"].resize(n_total, axis=0)
                for i, pb in enumerate(pose_blobs_raw):
                    grp["pose_blobs"][n_exist + i] = np.frombuffer(pb, dtype=np.uint8)

                grp["computed_at"].resize(n_total, axis=0)
                grp["computed_at"][n_exist:] = [c.encode("utf-8") for c in computed_ats]

                grp["source"].resize(n_total, axis=0)
                grp["source"][n_exist:] = [s.encode("utf-8") for s in sources]

                grp["top_n"].resize(n_total, axis=0)
                grp["top_n"][n_exist:] = top_ns
            else:
                chunk = min(n_new, 1024)
                grp.create_dataset(
                    "compound_hashes",
                    data=np.array([c.encode("utf-8") for c in compound_hashes], dtype=dt_str),
                    maxshape=(None,),
                    chunks=(chunk,),
                )
                grp.create_dataset(
                    "docking_scores",
                    data=scores,
                    maxshape=(None,),
                    chunks=(chunk,),
                )
                ds_poses = grp.create_dataset(
                    "pose_blobs",
                    shape=(n_new,),
                    dtype=vlen_bytes,
                    maxshape=(None,),
                    chunks=(chunk,),
                )
                for i, pb in enumerate(pose_blobs_raw):
                    ds_poses[i] = np.frombuffer(pb, dtype=np.uint8)

                grp.create_dataset(
                    "computed_at",
                    data=np.array([c.encode("utf-8") for c in computed_ats], dtype=dt_str),
                    maxshape=(None,),
                    chunks=(chunk,),
                )
                grp.create_dataset(
                    "source",
                    data=np.array([s.encode("utf-8") for s in sources], dtype=dt_str),
                    maxshape=(None,),
                    chunks=(chunk,),
                )
                grp.create_dataset(
                    "top_n",
                    data=top_ns,
                    maxshape=(None,),
                    chunks=(chunk,),
                )

            logger.info(f"write_bundle: {n_new} 件書き込み → protein_hash={protein_hash}")
            return n_new

    def read_bundle(self, protein_hash: str) -> Optional[List[Dict[str, Any]]]:
        """protein-bundle スキーマ (v3) から 1 タンパク質分の結果を読み込む。

        Returns list of dicts with keys: compound_hash, score, pose_blob.
        score == _FAILED_SCORE_SENTINEL (-999.0) は None に変換して返す。
        """
        if not self.hdf5_file_path.exists():
            return None

        group_path = f"/results/{protein_hash}"
        try:
            with h5py.File(self.hdf5_file_path, "r") as f:
                if group_path not in f or "compound_hashes" not in f[group_path]:
                    return None

                grp = f[group_path]
                n = grp["compound_hashes"].shape[0]
                raw_hashes = grp["compound_hashes"][:]
                compound_hashes_decoded = [
                    h.decode("utf-8") if isinstance(h, bytes) else str(h)
                    for h in raw_hashes
                ]
                docking_scores = grp["docking_scores"][:].tolist()
                pose_blobs_raw = [bytes(grp["pose_blobs"][i]) for i in range(n)]

                results = []
                for i, ch in enumerate(compound_hashes_decoded):
                    sc = docking_scores[i]
                    results.append({
                        "compound_hash": ch,
                        "score": None if sc == _FAILED_SCORE_SENTINEL else sc,
                        "pose_blob": pose_blobs_raw[i],
                    })
                return results

        except Exception as e:
            logger.error(f"read_bundle error (protein_hash={protein_hash}): {e}", exc_info=True)
            raise

    def _exists_bundle(self, protein_hash: str, compound_hash: str) -> bool:
        """protein-bundle スキーマ (v3) で指定ペアが存在するか確認。"""
        if not self.hdf5_file_path.exists():
            return False

        group_path = f"/results/{protein_hash}"
        try:
            with h5py.File(self.hdf5_file_path, "r") as f:
                if group_path not in f or "compound_hashes" not in f[group_path]:
                    return False
                hashes = {
                    h.decode("utf-8") if isinstance(h, bytes) else str(h)
                    for h in f[group_path]["compound_hashes"][:]
                }
                return compound_hash in hashes
        except Exception as e:
            logger.error(f"_exists_bundle error: {e}", exc_info=True)
            return False

    def get_all_keys_bundle(self) -> Set[Tuple[str, str]]:
        """protein-bundle スキーマ (v3) から全 (protein_hash, compound_hash) ペアを返す。"""
        keys: Set[Tuple[str, str]] = set()
        if not self.hdf5_file_path.exists():
            return keys

        try:
            with h5py.File(self.hdf5_file_path, "r") as f:
                if "results" not in f:
                    return keys
                for protein_hash in f["results"]:
                    grp = f["results"][protein_hash]
                    if "compound_hashes" not in grp:
                        continue
                    for raw_h in grp["compound_hashes"][:]:
                        ch = raw_h.decode("utf-8") if isinstance(raw_h, bytes) else str(raw_h)
                        keys.add((protein_hash, ch))
        except Exception as e:
            logger.error(f"get_all_keys_bundle error: {e}", exc_info=True)
            raise

        return keys
