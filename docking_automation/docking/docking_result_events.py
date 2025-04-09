"""
ドッキング結果に関連するドメインイベントを定義するモジュール。

このモジュールは、ドッキング結果の永続化や状態変更に関連するドメインイベントを提供します。
"""

from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, ClassVar, Dict, Optional

from docking_automation.docking.docking_result import DockingResult
from docking_automation.domain.domain_event import DomainEvent


class DockingResultEvent(DomainEvent):
    """
    ドッキング結果に関連するドメインイベントの基底クラス。

    ドッキング結果に関連する重要な出来事を表現するイベントオブジェクト。
    不変（イミュータブル）であり、過去に発生した事実を表現する。
    """

    event_type: ClassVar[str] = "DockingResultEvent"

    def __init__(self, result: DockingResult, occurred_on: Optional[datetime] = None):
        """
        DockingResultEventオブジェクトを初期化する。

        Args:
            result: ドッキング結果
            occurred_on: イベントが発生した日時（指定しない場合は現在時刻）
        """
        super().__init__(occurred_on=occurred_on or datetime.now())
        self._result = result

    @property
    def result(self) -> DockingResult:
        """
        ドッキング結果を取得する。

        Returns:
            ドッキング結果
        """
        return self._result

    def to_dict(self) -> Dict[str, Any]:
        """
        イベントを辞書形式に変換する。

        Returns:
            イベントの内容を表す辞書
        """
        result_dict = super().to_dict()
        result_dict["result_id"] = self.result.id
        result_dict["protein_id"] = self.result.protein_id
        result_dict["compound_set_id"] = self.result.compound_set_id
        result_dict["compound_index"] = self.result.compound_index
        result_dict["docking_score"] = self.result.docking_score
        return result_dict


class DockingResultCreatedEvent(DockingResultEvent):
    """
    ドッキング結果が作成されたことを表すイベント。
    """

    event_type: ClassVar[str] = "DockingResultCreatedEvent"


class DockingResultSavedEvent(DockingResultEvent):
    """
    ドッキング結果が保存されたことを表すイベント。
    """

    event_type: ClassVar[str] = "DockingResultSavedEvent"

    def __init__(self, result: DockingResult, storage_path: str, occurred_on: Optional[datetime] = None):
        """
        DockingResultSavedEventオブジェクトを初期化する。

        Args:
            result: ドッキング結果
            storage_path: 保存先のパス
            occurred_on: イベントが発生した日時（指定しない場合は現在時刻）
        """
        super().__init__(result=result, occurred_on=occurred_on)
        self._storage_path = storage_path

    @property
    def storage_path(self) -> str:
        """
        保存先のパスを取得する。

        Returns:
            保存先のパス
        """
        return self._storage_path

    def to_dict(self) -> Dict[str, Any]:
        """
        イベントを辞書形式に変換する。

        Returns:
            イベントの内容を表す辞書
        """
        result_dict = super().to_dict()
        result_dict["storage_path"] = self.storage_path
        return result_dict


class DockingResultLoadedEvent(DockingResultEvent):
    """
    ドッキング結果が読み込まれたことを表すイベント。
    """

    event_type: ClassVar[str] = "DockingResultLoadedEvent"

    def __init__(self, result: DockingResult, source_path: str, occurred_on: Optional[datetime] = None):
        """
        DockingResultLoadedEventオブジェクトを初期化する。

        Args:
            result: ドッキング結果
            source_path: 読み込み元のパス
            occurred_on: イベントが発生した日時（指定しない場合は現在時刻）
        """
        super().__init__(result=result, occurred_on=occurred_on)
        self._source_path = source_path

    @property
    def source_path(self) -> str:
        """
        読み込み元のパスを取得する。

        Returns:
            読み込み元のパス
        """
        return self._source_path

    def to_dict(self) -> Dict[str, Any]:
        """
        イベントを辞書形式に変換する。

        Returns:
            イベントの内容を表す辞書
        """
        result_dict = super().to_dict()
        result_dict["source_path"] = self.source_path
        return result_dict


class DockingResultUpdatedEvent(DockingResultEvent):
    """
    ドッキング結果が更新されたことを表すイベント。
    """

    event_type: ClassVar[str] = "DockingResultUpdatedEvent"

    def __init__(
        self, result: DockingResult, previous_version: int, new_version: int, occurred_on: Optional[datetime] = None
    ):
        """
        DockingResultUpdatedEventオブジェクトを初期化する。

        Args:
            result: ドッキング結果
            previous_version: 更新前のバージョン
            new_version: 更新後のバージョン
            occurred_on: イベントが発生した日時（指定しない場合は現在時刻）
        """
        super().__init__(result=result, occurred_on=occurred_on)
        self._previous_version = previous_version
        self._new_version = new_version

    @property
    def previous_version(self) -> int:
        """
        更新前のバージョンを取得する。

        Returns:
            更新前のバージョン
        """
        return self._previous_version

    @property
    def new_version(self) -> int:
        """
        更新後のバージョンを取得する。

        Returns:
            更新後のバージョン
        """
        return self._new_version

    def to_dict(self) -> Dict[str, Any]:
        """
        イベントを辞書形式に変換する。

        Returns:
            イベントの内容を表す辞書
        """
        result_dict = super().to_dict()
        result_dict["previous_version"] = self.previous_version
        result_dict["new_version"] = self.new_version
        return result_dict


class DockingResultDeletedEvent(DockingResultEvent):
    """
    ドッキング結果が削除されたことを表すイベント。
    """

    event_type: ClassVar[str] = "DockingResultDeletedEvent"
