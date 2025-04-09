"""
ハイブリッドドッキング結果リポジトリモジュール。

このモジュールは、SQLiteデータベースとファイルシステムを組み合わせたハイブリッドリポジトリを提供します。
メタデータはSQLiteデータベースに保存し、ポーズデータはファイルシステムに保存します。
"""

import json
import os
import shutil
import sqlite3
import tempfile
from contextlib import contextmanager
from datetime import datetime
from pathlib import Path
from threading import Lock
from typing import Any, Dict, Iterator, List, Optional, Set, Tuple, cast

import numpy as np
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


class HybridDockingResultRepository(DockingResultRepositoryABC):
    """
    ハイブリッドドッキング結果リポジトリ。

    SQLiteデータベースとファイルシステムを組み合わせたハイブリッドリポジトリです。
    メタデータはSQLiteデータベースに保存し、ポーズデータはファイルシステムに保存します。
    これにより、メタデータの効率的な検索と大容量ポーズデータの効率的な保存を両立します。
    """

    def __init__(self, database_path: Path, files_directory: Path):
        """
        HybridDockingResultRepositoryオブジェクトを初期化する。

        Args:
            database_path: データベースファイルのパス
            files_directory: ポーズデータファイルを保存するディレクトリ
        """
        self.database_path = database_path
        self.files_directory = files_directory
        self.locks_directory = files_directory.parent / "locks"

        # トランザクション管理用の変数
        self._transaction_lock = Lock()
        self._in_transaction = False
        self._connection: Optional[sqlite3.Connection] = None

        # ディレクトリ構造を作成
        self.files_directory.mkdir(parents=True, exist_ok=True)
        self.locks_directory.mkdir(parents=True, exist_ok=True)

        # データベースを初期化
        self._initialize_database()

        # トランザクション管理用の変数
        self._transaction_lock = Lock()
        self._in_transaction = False
        self._connection: Optional[sqlite3.Connection] = None

    def _initialize_database(self) -> None:
        """
        データベースを初期化する。
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # docking_resultsテーブルを作成
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS docking_results (
                    id TEXT PRIMARY KEY,
                    protein_id TEXT NOT NULL,
                    compound_set_id TEXT NOT NULL,
                    compound_index INTEGER NOT NULL,
                    docking_score REAL NOT NULL,
                    version INTEGER NOT NULL DEFAULT 1,
                    created_at TEXT NOT NULL,
                    updated_at TEXT NOT NULL,
                    UNIQUE(protein_id, compound_set_id, compound_index)
                )
            """
            )

            # result_filesテーブルを作成
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS result_files (
                    result_id TEXT PRIMARY KEY,
                    file_path TEXT NOT NULL,
                    file_type TEXT NOT NULL,
                    file_size INTEGER NOT NULL,
                    FOREIGN KEY(result_id) REFERENCES docking_results(id) ON DELETE CASCADE
                )
            """
            )

            # metadataテーブルを作成
            cursor.execute(
                """
                CREATE TABLE IF NOT EXISTS metadata (
                    result_id TEXT NOT NULL,
                    key TEXT NOT NULL,
                    value TEXT NOT NULL,
                    PRIMARY KEY(result_id, key),
                    FOREIGN KEY(result_id) REFERENCES docking_results(id) ON DELETE CASCADE
                )
            """
            )

            # インデックスを作成
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_docking_results_protein_id ON docking_results(protein_id)")
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_docking_results_compound_set_id ON docking_results(compound_set_id)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_docking_results_compound_index ON docking_results(compound_index)"
            )
            cursor.execute(
                "CREATE INDEX IF NOT EXISTS idx_docking_results_docking_score ON docking_results(docking_score)"
            )
            cursor.execute("CREATE INDEX IF NOT EXISTS idx_metadata_key ON metadata(key)")

            conn.commit()

    @contextmanager
    def _get_connection(self) -> Iterator[sqlite3.Connection]:
        """
        データベース接続を取得する。

        Yields:
            データベース接続
        """
        # トランザクション内の場合は既存の接続を使用
        if self._in_transaction and self._connection is not None:
            yield self._connection
            return

        # 新しい接続を作成
        conn = sqlite3.connect(str(self.database_path))

        # 外部キー制約を有効化
        conn.execute("PRAGMA foreign_keys = ON")

        # WAL（Write-Ahead Logging）モードを有効化
        conn.execute("PRAGMA journal_mode = WAL")

        # 行を辞書として取得できるようにする
        conn.row_factory = sqlite3.Row

        try:
            yield conn
            conn.commit()
        except Exception:
            conn.rollback()
            raise
        finally:
            if not self._in_transaction:
                conn.close()


def save(self, result: DockingResult) -> None:
    """
    ドッキング計算結果を保存する。

    Args:
        result: 保存する結果
    """
    # ポーズデータファイルを保存
    file_path = self._save_pose_file(result)

    # トランザクション内で実行されている場合は、既存の接続を使用
    if self._in_transaction and self._connection is not None:
        conn = self._connection
        cursor = conn.cursor()
        close_connection = False
    else:
        # 新しい接続を取得
        conn = self._get_connection()
        cursor = conn.cursor()
        close_connection = True

    try:
        # バージョンチェック（楽観的ロック）
        cursor.execute("SELECT version FROM docking_results WHERE id = ?", (result.id,))
        row = cursor.fetchone()

        if row is not None:
            current_version = row["version"]

            if current_version != result.version:
                raise ValueError(
                    f"バージョンの競合が発生しました: "
                    f"期待されるバージョン {result.version}, "
                    f"実際のバージョン {current_version}"
                )

            # 既存の結果を更新
            cursor.execute(
                """
                UPDATE docking_results
                SET protein_id = ?, compound_set_id = ?, compound_index = ?,
                    docking_score = ?, version = ?, updated_at = ?
                WHERE id = ?
                """,
                (
                    result.protein_id,
                    result.compound_set_id,
                    result.compound_index,
                    result.docking_score,
                    result.version + 1,
                    datetime.now().isoformat(),
                    result.id,
                ),
            )

            # 既存のメタデータを削除
            cursor.execute("DELETE FROM metadata WHERE result_id = ?", (result.id,))

            # 既存のファイル情報を更新
            cursor.execute(
                """
                UPDATE result_files
                SET file_path = ?, file_size = ?
                WHERE result_id = ?
                """,
                (
                    str(file_path),
                    os.path.getsize(file_path),
                    result.id,
                ),
            )

            # イベントを発行
            new_result = result.with_incremented_version()
            event = DockingResultUpdatedEvent(
                result=new_result,
                previous_version=result.version,
                new_version=new_result.version,
            )
            DomainEventPublisher.publish(event)
        else:
            try:
                # 新しい結果を挿入
                cursor.execute(
                    """
                    INSERT INTO docking_results
                    (id, protein_id, compound_set_id, compound_index, docking_score, version, created_at, updated_at)
                    VALUES (?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        result.id,
                        result.protein_id,
                        result.compound_set_id,
                        result.compound_index,
                        result.docking_score,
                        1,  # 新規作成時はバージョン1
                        datetime.now().isoformat(),
                        datetime.now().isoformat(),
                    ),
                )

                # ファイル情報を挿入
                cursor.execute(
                    """
                    INSERT INTO result_files
                    (result_id, file_path, file_type, file_size)
                    VALUES (?, ?, ?, ?)
                    """,
                    (
                        result.id,
                        str(file_path),
                        "sdf",
                        os.path.getsize(file_path),
                    ),
                )

                # メタデータを挿入
                for key, value in result.metadata.items():
                    # NumPy配列をリストに変換
                    serializable_value = self._make_json_serializable(value)
                    cursor.execute(
                        """
                        INSERT INTO metadata
                        (result_id, key, value)
                        VALUES (?, ?, ?)
                        """,
                        (
                            result.id,
                            key,
                            json.dumps(serializable_value),
                        ),
                    )

                # イベントを発行
                event = DockingResultSavedEvent(result=result, storage_path=str(file_path))
                DomainEventPublisher.publish(event)
            except sqlite3.IntegrityError as e:
                if "UNIQUE constraint failed" in str(e):
                    # 既存のレコードを更新
                    cursor.execute(
                        """
                        UPDATE docking_results
                        SET docking_score = ?, version = version + 1, updated_at = ?
                        WHERE protein_id = ? AND compound_set_id = ? AND compound_index = ?
                        """,
                        (
                            result.docking_score,
                            datetime.now().isoformat(),
                            result.protein_id,
                            result.compound_set_id,
                            result.compound_index,
                        ),
                    )

                    # 既存のメタデータを削除
                    cursor.execute("DELETE FROM metadata WHERE result_id = ?", (result.id,))

                    # メタデータを挿入
                    for key, value in result.metadata.items():
                        # NumPy配列をリストに変換
                        serializable_value = self._make_json_serializable(value)
                        cursor.execute(
                            """
                            INSERT INTO metadata
                            (result_id, key, value)
                            VALUES (?, ?, ?)
                            """,
                            (
                                result.id,
                                key,
                                json.dumps(serializable_value),
                            ),
                        )
                else:
                    raise
    finally:
        # トランザクション内で実行されていない場合は、接続を閉じる
        if close_connection:
            conn.close()

    def load(self, result_id: str) -> Optional[DockingResult]:
        """
        ドッキング計算結果を読み込む。

        Args:
            result_id: 結果のID

        Returns:
            読み込んだ結果、見つからなければNone
        """
        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 結果を取得
            cursor.execute(
                """
                SELECT r.*, f.file_path
                FROM docking_results r
                JOIN result_files f ON r.id = f.result_id
                WHERE r.id = ?
                """,
                (result_id,),
            )
            row = cursor.fetchone()

            if row is None:
                return None

            # メタデータを取得
            cursor.execute("SELECT key, value FROM metadata WHERE result_id = ?", (result_id,))
            metadata_rows = cursor.fetchall()

            metadata = {}
            for metadata_row in metadata_rows:
                metadata[metadata_row["key"]] = json.loads(metadata_row["value"])

            # DockingResultオブジェクトを作成
            result = DockingResult(
                result_path=Path(row["file_path"]),
                protein_id=row["protein_id"],
                compound_set_id=row["compound_set_id"],
                compound_index=row["compound_index"],
                docking_score=row["docking_score"],
                metadata=metadata,
                id=row["id"],
                version=row["version"],
            )

            # イベントを発行
            event = DockingResultLoadedEvent(result=result, source_path=row["file_path"])
            DomainEventPublisher.publish(event)

            return result

    def find_by_protein(self, protein_id: str) -> DockingResultCollection:
        """
        タンパク質IDに基づいて結果を検索する。

        Args:
            protein_id: タンパク質のID

        Returns:
            検索結果のコレクション
        """
        results: List[DockingResult] = []

        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 結果を取得
            cursor.execute(
                """
                SELECT r.*, f.file_path
                FROM docking_results r
                JOIN result_files f ON r.id = f.result_id
                WHERE r.protein_id = ?
                ORDER BY r.docking_score
                """,
                (protein_id,),
            )
            rows = cursor.fetchall()

            for row in rows:
                # メタデータを取得
                cursor.execute("SELECT key, value FROM metadata WHERE result_id = ?", (row["id"],))
                metadata_rows = cursor.fetchall()

                metadata = {}
                for metadata_row in metadata_rows:
                    metadata[metadata_row["key"]] = json.loads(metadata_row["value"])

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=Path(row["file_path"]),
                    protein_id=row["protein_id"],
                    compound_set_id=row["compound_set_id"],
                    compound_index=row["compound_index"],
                    docking_score=row["docking_score"],
                    metadata=metadata,
                    id=row["id"],
                    version=row["version"],
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

        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 結果を取得
            cursor.execute(
                """
                SELECT r.*, f.file_path
                FROM docking_results r
                JOIN result_files f ON r.id = f.result_id
                WHERE r.compound_set_id = ?
                ORDER BY r.docking_score
                """,
                (compound_set_id,),
            )
            rows = cursor.fetchall()

            for row in rows:
                # メタデータを取得
                cursor.execute("SELECT key, value FROM metadata WHERE result_id = ?", (row["id"],))
                metadata_rows = cursor.fetchall()

                metadata = {}
                for metadata_row in metadata_rows:
                    metadata[metadata_row["key"]] = json.loads(metadata_row["value"])

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=Path(row["file_path"]),
                    protein_id=row["protein_id"],
                    compound_set_id=row["compound_set_id"],
                    compound_index=row["compound_index"],
                    docking_score=row["docking_score"],
                    metadata=metadata,
                    id=row["id"],
                    version=row["version"],
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

        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 結果を取得
            cursor.execute(
                """
                SELECT r.*, f.file_path
                FROM docking_results r
                JOIN result_files f ON r.id = f.result_id
                WHERE r.compound_set_id = ? AND r.compound_index = ?
                ORDER BY r.docking_score
                """,
                (compound_set_id, compound_index),
            )
            rows = cursor.fetchall()

            for row in rows:
                # メタデータを取得
                cursor.execute("SELECT key, value FROM metadata WHERE result_id = ?", (row["id"],))
                metadata_rows = cursor.fetchall()

                metadata = {}
                for metadata_row in metadata_rows:
                    metadata[metadata_row["key"]] = json.loads(metadata_row["value"])

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=Path(row["file_path"]),
                    protein_id=row["protein_id"],
                    compound_set_id=row["compound_set_id"],
                    compound_index=row["compound_index"],
                    docking_score=row["docking_score"],
                    metadata=metadata,
                    id=row["id"],
                    version=row["version"],
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

        with self._get_connection() as conn:
            cursor = conn.cursor()

            # 結果を取得
            cursor.execute(
                """
                SELECT r.*, f.file_path
                FROM docking_results r
                JOIN result_files f ON r.id = f.result_id
                WHERE r.protein_id = ? AND r.compound_set_id = ? AND r.compound_index = ?
                ORDER BY r.docking_score
                """,
                (protein_id, compound_set_id, compound_index),
            )
            rows = cursor.fetchall()

            for row in rows:
                # メタデータを取得
                cursor.execute("SELECT key, value FROM metadata WHERE result_id = ?", (row["id"],))
                metadata_rows = cursor.fetchall()

                metadata = {}
                for metadata_row in metadata_rows:
                    metadata[metadata_row["key"]] = json.loads(metadata_row["value"])

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=Path(row["file_path"]),
                    protein_id=row["protein_id"],
                    compound_set_id=row["compound_set_id"],
                    compound_index=row["compound_index"],
                    docking_score=row["docking_score"],
                    metadata=metadata,
                    id=row["id"],
                    version=row["version"],
                )

                results.append(result)

        return DockingResultCollection(results)

    def save_collection(self, collection: DockingResultCollection) -> None:
        """
        ドッキング計算結果のコレクションを保存する。

        Args:
            collection: 保存するコレクション
        """
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
            self._connection = sqlite3.connect(str(self.database_path))

            # 外部キー制約を有効化
            self._connection.execute("PRAGMA foreign_keys = ON")

            # WAL（Write-Ahead Logging）モードを有効化
            self._connection.execute("PRAGMA journal_mode = WAL")

            # 行を辞書として取得できるようにする
            self._connection.row_factory = sqlite3.Row

            # トランザクションを開始
            self._connection.execute("BEGIN TRANSACTION")

    def commit_transaction(self) -> None:
        """
        トランザクションをコミットする。
        """
        with self._transaction_lock:
            if not self._in_transaction or self._connection is None:
                raise ValueError("トランザクションは開始されていません")

            # トランザクションをコミット
            self._connection.commit()

            # 接続を閉じる
            self._connection.close()
            self._connection = None

            # トランザクションをリセット
            self._in_transaction = False

    def rollback_transaction(self) -> None:
        """
        トランザクションをロールバックする。
        """
        with self._transaction_lock:
            if not self._in_transaction or self._connection is None:
                raise ValueError("トランザクションは開始されていません")

            # トランザクションをロールバック
            self._connection.rollback()

            # 接続を閉じる
            self._connection.close()
            self._connection = None

            # トランザクションをリセット
            self._in_transaction = False

    def _save_pose_file(self, result: DockingResult) -> Path:
        """
        ポーズデータファイルを保存する。

        Args:
            result: ドッキング結果

        Returns:
            保存されたファイルのパス
        """
        # ディレクトリ構造を作成
        result_dir = self.files_directory / result.protein_id / result.compound_set_id
        result_dir.mkdir(parents=True, exist_ok=True)

        # 結果ファイルのパス
        result_file = result_dir / f"{result.compound_index}.sdf"

        # ロックファイルのパス
        lock_file = self.locks_directory / f"{result.protein_id}_{result.compound_set_id}_{result.compound_index}.lock"

        # ファイルロックを取得
        with FileLock(str(lock_file)):
            # 一時ファイルに書き込み
            with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".tmp") as temp_file:
                # 結果ファイルをコピー
                shutil.copy2(result.result_path, temp_file.name)

            # 一時ファイルを本来のファイルに移動（アトミック操作）
            shutil.move(temp_file.name, result_file)

        return result_file

    def _make_json_serializable(self, obj: Any) -> Any:
        """
        オブジェクトをJSON化可能な形式に変換する。

        Args:
            obj: 変換対象のオブジェクト

        Returns:
            JSON化可能な形式に変換されたオブジェクト
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        elif isinstance(obj, dict):
            return {k: self._make_json_serializable(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [self._make_json_serializable(item) for item in obj]
        elif isinstance(obj, tuple):
            return tuple(self._make_json_serializable(item) for item in obj)
        else:
            return obj

    @classmethod
    def create(
        cls, database_path: Optional[Path] = None, files_directory: Optional[Path] = None
    ) -> "HybridDockingResultRepository":
        """
        HybridDockingResultRepositoryオブジェクトを作成するファクトリメソッド。

        Args:
            database_path: データベースファイルのパス（指定しない場合はデフォルトのパスを使用）
            files_directory: ポーズデータファイルを保存するディレクトリ（指定しない場合はデフォルトのパスを使用）

        Returns:
            作成されたHybridDockingResultRepositoryオブジェクト
        """
        if database_path is None:
            # デフォルトのパスを使用
            database_path = Path.home() / ".docking_automation" / "repository" / "docking_results.db"

        if files_directory is None:
            # デフォルトのパスを使用
            files_directory = Path.home() / ".docking_automation" / "repository" / "files"

        # データベースディレクトリを作成
        database_path.parent.mkdir(parents=True, exist_ok=True)

        return cls(database_path=database_path, files_directory=files_directory)
