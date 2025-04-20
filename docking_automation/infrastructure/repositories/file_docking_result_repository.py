"""
ファイルベースのドッキング結果リポジトリモジュール。

このモジュールは、ファイルシステムを使用してドッキング結果を永続化するリポジトリを提供します。
"""

import json
import os
import shutil
import tempfile
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from threading import Lock
from typing import Any, Dict, List, Optional, Set, Tuple, cast

from filelock import FileLock

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.docking.docking_result_events import (
    DockingResultCreatedEvent,
    DockingResultDeletedEvent,
    DockingResultLoadedEvent,
    DockingResultSavedEvent,
    DockingResultUpdatedEvent,
)
from docking_automation.domain.domain_event_publisher import DomainEventPublisher
from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepositoryABC,
)


class FileDockingResultRepository(DockingResultRepositoryABC):
    """
    ファイルベースのドッキング結果リポジトリ。

    ファイルシステムを使用してドッキング結果を永続化します。
    並列アクセスの競合を防ぐために、ファイルロック機構と一時ファイル + アトミック移動を使用します。
    """

    def __init__(self, base_directory: Path):
        """
        FileDockingResultRepositoryオブジェクトを初期化する。

        Args:
            base_directory: リポジトリのベースディレクトリ
        """
        self.base_directory = base_directory
        self.results_directory = base_directory / "results"
        self.metadata_directory = base_directory / "metadata"
        self.index_directory = base_directory / "index"
        self.locks_directory = base_directory / "locks"

        # ディレクトリ構造を作成
        self.results_directory.mkdir(parents=True, exist_ok=True)
        self.metadata_directory.mkdir(parents=True, exist_ok=True)
        self.index_directory.mkdir(parents=True, exist_ok=True)
        self.locks_directory.mkdir(parents=True, exist_ok=True)

        # トランザクション管理用の変数
        self._transaction_lock = Lock()
        self._in_transaction = False
        self._transaction_operations: List[Tuple[str, Any]] = []

        # スレッドプール
        self._executor = ThreadPoolExecutor(max_workers=10)

    def save(self, result: DockingResult) -> None:
        """
        ドッキング計算結果を保存する。

        Args:
            result: 保存する結果
        """
        # トランザクション内の場合は操作をキューに追加
        if self._in_transaction:
            self._transaction_operations.append(("save", result))
            return

        # 結果ディレクトリを作成
        result_dir = self._get_result_directory(result)
        result_dir.mkdir(parents=True, exist_ok=True)

        # メタデータディレクトリを作成
        metadata_dir = self._get_metadata_directory(result)
        metadata_dir.mkdir(parents=True, exist_ok=True)

        # 結果ファイルのパス
        result_file = self._get_result_file_path(result)
        metadata_file = self._get_metadata_file_path(result)

        # ロックファイルのパス
        lock_file = self._get_lock_file_path(result)

        # ファイルロックを取得
        with FileLock(str(lock_file)):
            # バージョンチェック（楽観的ロック）
            if result_file.exists() and metadata_file.exists():
                current_metadata = self._load_metadata(metadata_file)
                current_version = current_metadata.get("version", 1)

                if current_version != result.version:
                    raise ValueError(
                        f"バージョンの競合が発生しました: "
                        f"期待されるバージョン {result.version}, "
                        f"実際のバージョン {current_version}"
                    )

            # 一時ファイルに書き込み
            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_result_file:
                # 結果ファイルをコピー
                shutil.copy2(result.result_path, temp_result_file.name)

            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_metadata_file:
                # メタデータを作成
                metadata = {
                    "id": result.id,
                    "protein_id": result.protein_id,
                    "compound_set_id": result.compound_set_id,
                    "compound_index": result.compound_index,
                    "docking_score": result.docking_score,
                    "protein_content_hash": result.protein_content_hash,
                    "compoundset_content_hash": result.compoundset_content_hash,
                    "metadata": result.metadata,
                    "version": result.version + 1,  # バージョンをインクリメント
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                }

                # メタデータをJSONとして書き込み
                json.dump(metadata, temp_metadata_file, indent=2)

            # 一時ファイルを本来のファイルに移動（アトミック操作）
            shutil.move(temp_result_file.name, result_file)
            shutil.move(temp_metadata_file.name, metadata_file)

            # インデックスを更新
            self._update_index(result)

            # イベントを発行
            new_result = result.with_incremented_version()
            event = DockingResultSavedEvent(result=new_result, storage_path=str(result_file))
            DomainEventPublisher.publish(event)

    def load(self, result_id: str) -> Optional[DockingResult]:
        """
        ドッキング計算結果を読み込む。

        Args:
            result_id: 結果のID

        Returns:
            読み込んだ結果、見つからなければNone
        """
        # IDからパスを特定
        for protein_dir in self.metadata_directory.iterdir():
            if not protein_dir.is_dir():
                continue

            for compound_set_dir in protein_dir.iterdir():
                if not compound_set_dir.is_dir():
                    continue

                for metadata_file in compound_set_dir.iterdir():
                    if not metadata_file.is_file() or metadata_file.suffix != ".json":
                        continue

                    # メタデータを読み込み
                    metadata = self._load_metadata(metadata_file)

                    if metadata.get("id") == result_id:
                        # 対応する結果ファイルのパス
                        result_file = self._get_result_file_path_from_metadata(metadata)

                        if not result_file.exists():
                            return None

                        # DockingResultオブジェクトを作成
                        result = DockingResult(
                            result_path=result_file,
                            protein_id=metadata["protein_id"],
                            compound_set_id=metadata["compound_set_id"],
                            compound_index=metadata["compound_index"],
                            docking_score=metadata["docking_score"],
                            protein_content_hash=metadata.get("protein_content_hash", "default_protein_hash"),
                            compoundset_content_hash=metadata.get("compoundset_content_hash", "default_compound_hash"),
                            metadata=metadata.get("metadata", {}),
                            id=metadata["id"],
                            version=metadata.get("version", 1),
                        )

                        # イベントを発行
                        event = DockingResultLoadedEvent(result=result, source_path=str(result_file))
                        DomainEventPublisher.publish(event)

                        return result

        return None

    def find_by_protein(self, protein_id: str) -> DockingResultCollection:
        """
        タンパク質IDに基づいて結果を検索する。

        Args:
            protein_id: タンパク質のID

        Returns:
            検索結果のコレクション
        """
        results: List[DockingResult] = []

        # タンパク質IDに対応するディレクトリを検索
        protein_dir = self.metadata_directory / protein_id

        if not protein_dir.exists() or not protein_dir.is_dir():
            return DockingResultCollection(results)

        # すべての化合物セットディレクトリを検索
        for compound_set_dir in protein_dir.iterdir():
            if not compound_set_dir.is_dir():
                continue

            # すべてのメタデータファイルを検索
            for metadata_file in compound_set_dir.iterdir():
                if not metadata_file.is_file() or metadata_file.suffix != ".json":
                    continue

                # メタデータを読み込み
                metadata = self._load_metadata(metadata_file)

                # 対応する結果ファイルのパス
                result_file = self._get_result_file_path_from_metadata(metadata)

                if not result_file.exists():
                    continue

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=result_file,
                    protein_id=metadata["protein_id"],
                    compound_set_id=metadata["compound_set_id"],
                    compound_index=metadata["compound_index"],
                    docking_score=metadata["docking_score"],
                    protein_content_hash=metadata.get("protein_content_hash", "default_protein_hash"),
                    compoundset_content_hash=metadata.get("compoundset_content_hash", "default_compound_hash"),
                    metadata=metadata.get("metadata", {}),
                    id=metadata["id"],
                    version=metadata.get("version", 1),
                )

                results.append(result)

        return DockingResultCollection(results)

    def find_by_compound_set(self, compound_set_id: str) -> DockingResultCollection:
        """
        化合物セットIDに基づいて結果を検索する。

        Args:
            compound_set_id: 化合物セットのID

        Returns:
            検索結果のコレクション
        """
        results: List[DockingResult] = []

        # すべてのタンパク質ディレクトリを検索
        for protein_dir in self.metadata_directory.iterdir():
            if not protein_dir.is_dir():
                continue

            # 化合物セットIDに対応するディレクトリを検索
            compound_set_dir = protein_dir / compound_set_id

            if not compound_set_dir.exists() or not compound_set_dir.is_dir():
                continue

            # すべてのメタデータファイルを検索
            for metadata_file in compound_set_dir.iterdir():
                if not metadata_file.is_file() or metadata_file.suffix != ".json":
                    continue

                # メタデータを読み込み
                metadata = self._load_metadata(metadata_file)

                # 対応する結果ファイルのパス
                result_file = self._get_result_file_path_from_metadata(metadata)

                if not result_file.exists():
                    continue

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=result_file,
                    protein_id=metadata["protein_id"],
                    compound_set_id=metadata["compound_set_id"],
                    compound_index=metadata["compound_index"],
                    docking_score=metadata["docking_score"],
                    protein_content_hash=metadata.get("protein_content_hash", "default_protein_hash"),
                    compoundset_content_hash=metadata.get("compoundset_content_hash", "default_compound_hash"),
                    metadata=metadata.get("metadata", {}),
                    id=metadata["id"],
                    version=metadata.get("version", 1),
                )

                results.append(result)

        return DockingResultCollection(results)

    def find_by_compound(self, compound_set_id: str, compound_index: int) -> DockingResultCollection:
        """
        化合物（化合物セットIDと化合物インデックスの組）に基づいて結果を検索する。

        Args:
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス

        Returns:
            検索結果のコレクション
        """
        results: List[DockingResult] = []

        # すべてのタンパク質ディレクトリを検索
        for protein_dir in self.metadata_directory.iterdir():
            if not protein_dir.is_dir():
                continue

            # 化合物セットIDに対応するディレクトリを検索
            compound_set_dir = protein_dir / compound_set_id

            if not compound_set_dir.exists() or not compound_set_dir.is_dir():
                continue

            # 化合物インデックスに対応するメタデータファイルを検索
            metadata_file = compound_set_dir / f"{compound_index}.json"

            if not metadata_file.exists() or not metadata_file.is_file():
                continue

            # メタデータを読み込み
            metadata = self._load_metadata(metadata_file)

            # 対応する結果ファイルのパス
            result_file = self._get_result_file_path_from_metadata(metadata)

            if not result_file.exists():
                continue

            # DockingResultオブジェクトを作成
            result = DockingResult(
                result_path=result_file,
                protein_id=metadata["protein_id"],
                compound_set_id=metadata["compound_set_id"],
                compound_index=metadata["compound_index"],
                docking_score=metadata["docking_score"],
                protein_content_hash=metadata.get("protein_content_hash", "default_protein_hash"),
                compoundset_content_hash=metadata.get("compoundset_content_hash", "default_compound_hash"),
                metadata=metadata.get("metadata", {}),
                id=metadata["id"],
                version=metadata.get("version", 1),
            )

            results.append(result)

        return DockingResultCollection(results)

    def find_by_protein_and_compound(
        self, protein_id: str, compound_set_id: str, compound_index: int
    ) -> DockingResultCollection:
        """
        タンパク質IDと化合物（化合物セットIDと化合物インデックスの組）に基づいて結果を検索する。

        Args:
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス

        Returns:
            検索結果のコレクション
        """
        results: List[DockingResult] = []

        # タンパク質IDに対応するディレクトリを検索
        protein_dir = self.metadata_directory / protein_id

        if not protein_dir.exists() or not protein_dir.is_dir():
            return DockingResultCollection(results)

        # 化合物セットIDに対応するディレクトリを検索
        compound_set_dir = protein_dir / compound_set_id

        if not compound_set_dir.exists() or not compound_set_dir.is_dir():
            return DockingResultCollection(results)

        # 化合物インデックスに対応するメタデータファイルを検索
        metadata_file = compound_set_dir / f"{compound_index}.json"

        if not metadata_file.exists() or not metadata_file.is_file():
            return DockingResultCollection(results)

        # メタデータを読み込み
        metadata = self._load_metadata(metadata_file)

        # 対応する結果ファイルのパス
        result_file = self._get_result_file_path_from_metadata(metadata)

        if not result_file.exists():
            return DockingResultCollection(results)

        # DockingResultオブジェクトを作成
        result = DockingResult(
            result_path=result_file,
            protein_id=metadata["protein_id"],
            compound_set_id=metadata["compound_set_id"],
            compound_index=metadata["compound_index"],
            docking_score=metadata["docking_score"],
            protein_content_hash=metadata.get("protein_content_hash", "default_protein_hash"),
            compoundset_content_hash=metadata.get("compoundset_content_hash", "default_compound_hash"),
            metadata=metadata.get("metadata", {}),
            id=metadata["id"],
            version=metadata.get("version", 1),
        )

        results.append(result)

        return DockingResultCollection(results)

    def save_collection(self, collection: DockingResultCollection) -> None:
        """
        ドッキング計算結果のコレクションを保存する。

        Args:
            collection: 保存するコレクション
        """
        # トランザクション内の場合は操作をキューに追加
        if self._in_transaction:
            self._transaction_operations.append(("save_collection", collection))
            return

        # トランザクションを開始
        self.begin_transaction()

        try:
            # コレクション内のすべての結果を保存
            for result in collection:
                self.save(result)

            # トランザクションをコミット
            self.commit_transaction()
        except Exception as e:
            # エラーが発生した場合はロールバック
            self.rollback_transaction()
            raise e

    def begin_transaction(self) -> None:
        """
        トランザクションを開始する。
        """
        with self._transaction_lock:
            if self._in_transaction:
                raise ValueError("トランザクションは既に開始されています")

            self._in_transaction = True
            self._transaction_operations = []

    def commit_transaction(self) -> None:
        """
        トランザクションをコミットする。
        """
        with self._transaction_lock:
            if not self._in_transaction:
                raise ValueError("トランザクションは開始されていません")

            # トランザクション内の操作を実行
            for operation, args in self._transaction_operations:
                if operation == "save":
                    self._save_without_transaction(cast(DockingResult, args))
                elif operation == "save_collection":
                    for result in cast(DockingResultCollection, args):
                        self._save_without_transaction(result)

            # トランザクションをリセット
            self._in_transaction = False
            self._transaction_operations = []

    def rollback_transaction(self) -> None:
        """
        トランザクションをロールバックする。
        """
        with self._transaction_lock:
            if not self._in_transaction:
                raise ValueError("トランザクションは開始されていません")

            # トランザクションをリセット
            self._in_transaction = False
            self._transaction_operations = []

    def _save_without_transaction(self, result: DockingResult) -> None:
        """
        トランザクションなしでドッキング計算結果を保存する。

        Args:
            result: 保存する結果
        """
        # 結果ディレクトリを作成
        result_dir = self._get_result_directory(result)
        result_dir.mkdir(parents=True, exist_ok=True)

        # メタデータディレクトリを作成
        metadata_dir = self._get_metadata_directory(result)
        metadata_dir.mkdir(parents=True, exist_ok=True)

        # 結果ファイルのパス
        result_file = self._get_result_file_path(result)
        metadata_file = self._get_metadata_file_path(result)

        # ロックファイルのパス
        lock_file = self._get_lock_file_path(result)

        # ファイルロックを取得
        with FileLock(str(lock_file)):
            # バージョンチェック（楽観的ロック）
            if result_file.exists() and metadata_file.exists():
                current_metadata = self._load_metadata(metadata_file)
                current_version = current_metadata.get("version", 1)

                if current_version != result.version:
                    raise ValueError(
                        f"バージョンの競合が発生しました: "
                        f"期待されるバージョン {result.version}, "
                        f"実際のバージョン {current_version}"
                    )

            # 一時ファイルに書き込み
            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_result_file:
                # 結果ファイルをコピー
                shutil.copy2(result.result_path, temp_result_file.name)

            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_metadata_file:
                # メタデータを作成
                metadata = {
                    "id": result.id,
                    "protein_id": result.protein_id,
                    "compound_set_id": result.compound_set_id,
                    "compound_index": result.compound_index,
                    "docking_score": result.docking_score,
                    "protein_content_hash": result.protein_content_hash,
                    "compoundset_content_hash": result.compoundset_content_hash,
                    "metadata": result.metadata,
                    "version": result.version + 1,  # バージョンをインクリメント
                    "created_at": datetime.now().isoformat(),
                    "updated_at": datetime.now().isoformat(),
                }

                # メタデータをJSONとして書き込み
                json.dump(metadata, temp_metadata_file, indent=2)

            # 一時ファイルを本来のファイルに移動（アトミック操作）
            shutil.move(temp_result_file.name, result_file)
            shutil.move(temp_metadata_file.name, metadata_file)

            # インデックスを更新
            self._update_index(result)

            # イベントを発行
            new_result = result.with_incremented_version()
            event = DockingResultSavedEvent(result=new_result, storage_path=str(result_file))
            DomainEventPublisher.publish(event)

    def _get_result_directory(self, result: DockingResult) -> Path:
        """
        結果ディレクトリのパスを取得する。

        Args:
            result: ドッキング結果

        Returns:
            結果ディレクトリのパス
        """
        return self.results_directory / result.protein_id / result.compound_set_id

    def _get_metadata_directory(self, result: DockingResult) -> Path:
        """
        メタデータディレクトリのパスを取得する。

        Args:
            result: ドッキング結果

        Returns:
            メタデータディレクトリのパス
        """
        return self.metadata_directory / result.protein_id / result.compound_set_id

    def _get_result_file_path(self, result: DockingResult) -> Path:
        """
        結果ファイルのパスを取得する。

        Args:
            result: ドッキング結果

        Returns:
            結果ファイルのパス
        """
        return self._get_result_directory(result) / f"{result.compound_index}.sdf"

    def _get_metadata_file_path(self, result: DockingResult) -> Path:
        """
        メタデータファイルのパスを取得する。

        Args:
            result: ドッキング結果

        Returns:
            メタデータファイルのパス
        """
        return self._get_metadata_directory(result) / f"{result.compound_index}.json"

    def _get_lock_file_path(self, result: DockingResult) -> Path:
        """
        ロックファイルのパスを取得する。

        Args:
            result: ドッキング結果

        Returns:
            ロックファイルのパス
        """
        return self.locks_directory / f"{result.protein_id}_{result.compound_set_id}_{result.compound_index}.lock"

    def _get_result_file_path_from_metadata(self, metadata: Dict[str, Any]) -> Path:
        """
        メタデータから結果ファイルのパスを取得する。

        Args:
            metadata: メタデータ

        Returns:
            結果ファイルのパス
        """
        return (
            self.results_directory
            / metadata["protein_id"]
            / metadata["compound_set_id"]
            / f"{metadata['compound_index']}.sdf"
        )

    def _load_metadata(self, metadata_file: Path) -> Dict[str, Any]:
        """
        メタデータファイルを読み込む。

        Args:
            metadata_file: メタデータファイルのパス

        Returns:
            メタデータ
        """
        with open(metadata_file, "r") as f:
            return json.load(f)

    def _update_index(self, result: DockingResult) -> None:
        """
        インデックスを更新する。

        Args:
            result: ドッキング結果
        """
        # タンパク質IDによるインデックスを更新
        protein_index_file = self.index_directory / f"{result.protein_id}.idx"
        protein_index_lock = self.locks_directory / f"{result.protein_id}.idx.lock"

        with FileLock(str(protein_index_lock)):
            protein_index: Dict[str, List[str]] = {}

            if protein_index_file.exists():
                with open(protein_index_file, "r") as f:
                    protein_index = json.load(f)

            compound_key = f"{result.compound_set_id}_{result.compound_index}"

            if compound_key not in protein_index.get("compounds", []):
                if "compounds" not in protein_index:
                    protein_index["compounds"] = []

                protein_index["compounds"].append(compound_key)

            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_file:
                json.dump(protein_index, temp_file, indent=2)

            shutil.move(temp_file.name, protein_index_file)

        # 化合物セットIDによるインデックスを更新
        compound_set_index_file = self.index_directory / f"{result.compound_set_id}.idx"
        compound_set_index_lock = self.locks_directory / f"{result.compound_set_id}.idx.lock"

        with FileLock(str(compound_set_index_lock)):
            compound_set_index: Dict[str, List[str]] = {}

            if compound_set_index_file.exists():
                with open(compound_set_index_file, "r") as f:
                    compound_set_index = json.load(f)

            protein_key = f"{result.protein_id}_{result.compound_index}"

            if protein_key not in compound_set_index.get("proteins", []):
                if "proteins" not in compound_set_index:
                    compound_set_index["proteins"] = []

                compound_set_index["proteins"].append(protein_key)

            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_file:
                json.dump(compound_set_index, temp_file, indent=2)

            shutil.move(temp_file.name, compound_set_index_file)

    @classmethod
    def create(cls, base_directory: Optional[Path] = None) -> "FileDockingResultRepository":
        """
        FileDockingResultRepositoryオブジェクトを作成するファクトリメソッド。

        Args:
            base_directory: リポジトリのベースディレクトリ（指定しない場合はデフォルトのパスを使用）

        Returns:
            作成されたFileDockingResultRepositoryオブジェクト
        """
        if base_directory is None:
            # デフォルトのパスを使用
            base_directory = Path.home() / ".docking_automation" / "repository"

        return cls(base_directory=base_directory)
