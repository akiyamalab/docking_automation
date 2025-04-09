from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Set

from docking_automation.domain.domain_event import DomainEvent

# - 属性値に基づくハッシュ値計算を実装する
# - 値の変更が必要な場合は新しいインスタンスを返すメソッドを提供する
# - DockingResultCollectionのために、Sortableでなければならない。ドッキングスコアでソートするようにする。
# - result_pathに存在する分子構造は、ただ1つのポーズを持つことを保証する。


# エンティティ
class DockingResult:
    """
    ドッキング計算の結果を保持するクラス。
    1つのタンパク質と1つの化合物のドッキング計算結果が保持される。

    TODO: [DDD] エンティティとしての振る舞いを強化する
    - 一意の識別子（ID）を追加する
    - 不変属性と可変属性を明確に分離する
    - ライフサイクルを明確に定義する
    """

    def __init__(
        self,
        result_path: Path,
        protein_id: str,
        compound_set_id: str,
        compound_index: int,
        docking_score: float,
        metadata: Optional[Dict[str, Any]] = None,
        id: Optional[str] = None,
        version: int = 1,
    ):
        """
        DockingResultオブジェクトを初期化する。

        Args:
            result_path: 結果 SDF ファイルのパス
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            docking_score: ドッキングスコア
            metadata: メタデータ
            id: 結果のID（指定しない場合は自動生成）
            version: 結果のバージョン（楽観的ロックのために使用）
        """
        # 不変条件のバリデーション
        if not isinstance(result_path, Path):
            raise ValueError("result_pathはPathオブジェクトである必要があります")
        if not protein_id:
            raise ValueError("protein_idは空であってはなりません")
        if not compound_set_id:
            raise ValueError("compound_set_idは空であってはなりません")
        if compound_index < 0:
            raise ValueError("compound_indexは0以上である必要があります")
        if version < 1:
            raise ValueError("versionは1以上である必要があります")

        self.result_path = result_path
        self.protein_id = protein_id
        self.compound_set_id = compound_set_id
        self.compound_index = compound_index
        self.docking_score = docking_score
        self.metadata = metadata or {}
        self.id = id or f"{protein_id}_{compound_set_id}_{compound_index}_{datetime.now().strftime('%Y%m%d%H%M%S')}"
        self.version = version
        self._domain_events: Set[DomainEvent] = set()

    def __eq__(self, other: object) -> bool:
        """
        等価性比較を行う。DockingResultの場合はIDで比較する。

        Args:
            other: 比較対象のオブジェクト

        Returns:
            等価であればTrue、そうでなければFalse
        """
        if not isinstance(other, DockingResult):
            return False
        return self.id == other.id

    def __hash__(self) -> int:
        """
        ハッシュ値を計算する。DockingResultの場合はIDを使用する。

        Returns:
            ハッシュ値
        """
        return hash(self.id)

    def get_pose(self) -> Path:
        """
        ドッキングポーズのファイルパスを取得する。

        Returns:
            ポーズファイルのパス
        """
        return self.result_path

    # TODO: [DDD] ドメインロジックを追加する
    def get_score(self) -> float:
        """
        ドッキングスコアを取得する。

        Returns:
            ドッキングスコア
        """
        return self.docking_score

    # TODO: [DDD] ドメインロジックを追加する
    def get_metadata_value(self, key: str, default: Any = None) -> Any:
        """
        メタデータから指定されたキーの値を取得する。

        Args:
            key: 取得するメタデータのキー
            default: キーが存在しない場合のデフォルト値

        Returns:
            メタデータの値
        """
        return self.metadata.get(key, default)

    def get_original_index(self, compound_sets: List[Any]) -> int:
        """
        元の化合物インデックスを計算する。

        分割された化合物セットの場合、元のインデックスを計算して返す。

        Args:
            compound_sets: 化合物セットのリスト

        Returns:
            元の化合物インデックス
        """
        original_index = self.compound_index

        # 化合物セットのIDが一致するものを探す
        for compound_set in compound_sets:
            if compound_set.id == self.compound_set_id:
                # インデックス範囲を取得
                properties = compound_set.get_properties()
                if "index_range" in properties:
                    index_range = properties["index_range"]
                    # 元のインデックスを計算
                    original_index = index_range["start"] + self.compound_index
                    break

        return original_index

    def format_result(self, rank: Optional[int] = None, original_index: Optional[int] = None) -> str:
        """
        結果を整形して文字列として返す。

        Args:
            rank: 結果のランク（指定しない場合は表示しない）
            original_index: 元の化合物インデックス（指定しない場合はcompound_indexを使用）

        Returns:
            整形された結果文字列
        """
        # 元のインデックスが指定されていない場合は、compound_indexを使用
        if original_index is None:
            original_index = self.compound_index

        # ランクが指定されている場合は、ランクを表示
        rank_str = f"{rank}. " if rank is not None else ""

        # 結果を整形
        result_str = f"{rank_str}Score: {self.docking_score}, Compound: actives_subset_{original_index}"

        return result_str

    def get_compound_info(self, original_index: Optional[int] = None) -> str:
        """
        化合物の詳細情報を文字列として返す。

        Args:
            original_index: 元の化合物インデックス（指定しない場合はcompound_indexを使用）

        Returns:
            化合物の詳細情報
        """
        # 元のインデックスが指定されていない場合は、compound_indexを使用
        if original_index is None:
            original_index = self.compound_index

        # メタデータに化合物名がある場合は、それを返す
        if "compound_name" in self.metadata:
            return f"化合物名: {self.metadata['compound_name']}"
        else:
            # メタデータに化合物名がない場合は、化合物IDを返す
            return f"化合物ID: actives_subset_{original_index}"

    # TODO: [DDD] 値オブジェクトの不変性を保持するメソッドを実装する
    def with_metadata(self, key: str, value: Any) -> "DockingResult":
        """
        新しいメタデータを持つ新しいDockingResultインスタンスを作成する。
        値オブジェクトは不変なので、既存のインスタンスを変更する代わりに新しいインスタンスを返す。

        Args:
            key: 追加するメタデータのキー
            value: 追加するメタデータの値

        Returns:
            新しいメタデータを持つ新しいDockingResultインスタンス
        """
        new_metadata = self.metadata.copy()
        new_metadata[key] = value

        return DockingResult(
            result_path=self.result_path,
            protein_id=self.protein_id,
            compound_set_id=self.compound_set_id,
            compound_index=self.compound_index,
            docking_score=self.docking_score,
            metadata=new_metadata,
        )

    @classmethod
    def create(
        cls,
        result_path: Path,
        protein_id: str,
        compound_set_id: str,
        compound_index: int,
        docking_score: float,
        metadata: Optional[Dict[str, Any]] = None,
        id: Optional[str] = None,
    ) -> "DockingResult":
        """
        DockingResultオブジェクトを作成するファクトリメソッド。
        値オブジェクトの生成を一元管理し、バリデーションを行う。

        Args:
            result_path: 結果 SDF ファイルのパス
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            docking_score: ドッキングスコア
            metadata: メタデータ
            id: 結果のID（指定しない場合は自動生成）

        Returns:
            作成されたDockingResultオブジェクト
        """
        # 新しいDockingResultオブジェクトを作成
        result = cls(
            result_path=result_path,
            protein_id=protein_id,
            compound_set_id=compound_set_id,
            compound_index=compound_index,
            docking_score=docking_score,
            metadata=metadata,
            id=id,
            version=1,  # 新規作成時はバージョン1
        )

        # DockingResultCreatedEventを登録
        from docking_automation.docking.docking_result_events import (
            DockingResultCreatedEvent,
        )

        result.register_domain_event(DockingResultCreatedEvent(result=result))

        return result

    def register_domain_event(self, event: DomainEvent) -> None:
        """
        ドメインイベントを登録する。

        Args:
            event: 登録するドメインイベント
        """
        self._domain_events.add(event)

    def clear_domain_events(self) -> None:
        """
        登録されたドメインイベントをクリアする。
        """
        self._domain_events.clear()

    def get_domain_events(self) -> Set[DomainEvent]:
        """
        登録されたドメインイベントを取得する。

        Returns:
            登録されたドメインイベントのセット
        """
        return self._domain_events.copy()

    def with_domain_event(self, event: DomainEvent) -> "DockingResult":
        """
        ドメインイベントを含む新しいDockingResultインスタンスを作成する。
        値オブジェクトは不変なので、イベントを登録した新しいインスタンスを返す。

        Args:
            event: 登録するドメインイベント

        Returns:
            ドメインイベントを含む新しいDockingResultインスタンス
        """
        # 新しいインスタンスを作成
        new_result = DockingResult(
            result_path=self.result_path,
            protein_id=self.protein_id,
            compound_set_id=self.compound_set_id,
            compound_index=self.compound_index,
            docking_score=self.docking_score,
            metadata=self.metadata.copy(),
            id=self.id,
            version=self.version,
        )

        # 既存のドメインイベントをコピー
        for existing_event in self._domain_events:
            new_result.register_domain_event(existing_event)

        # 新しいイベントを登録
        new_result.register_domain_event(event)

        return new_result

    def with_incremented_version(self) -> "DockingResult":
        """
        バージョンをインクリメントした新しいDockingResultインスタンスを作成する。
        楽観的ロックのために使用する。

        Returns:
            バージョンをインクリメントした新しいDockingResultインスタンス
        """
        return DockingResult(
            result_path=self.result_path,
            protein_id=self.protein_id,
            compound_set_id=self.compound_set_id,
            compound_index=self.compound_index,
            docking_score=self.docking_score,
            metadata=self.metadata.copy(),
            id=self.id,
            version=self.version + 1,
        )
