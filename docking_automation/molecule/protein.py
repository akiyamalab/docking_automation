import uuid
from pathlib import Path
from typing import Any, Dict, List, Optional

from docking_automation.domain.domain_event import DomainEvent
from docking_automation.molecule.protein_events import (
    ProteinPathUpdated,
    ProteinRegistered,
)

# TODO: [P1] [DDD] リポジトリパターンを導入し、Proteinエンティティの永続化を担当するリポジトリを実装する


class Protein:
    """
    タンパク質を表すクラス。
    現時点では .pdb 形式をサポートしている。

    IDを持つ不変エンティティとして実装されており、初期化後に属性を変更することはできない。
    状態変化が必要な場合は、新しいインスタンスを作成する。
    """

    def __init__(self, path: str | Path, id: Optional[str] = None):
        """
        Proteinオブジェクトを初期化する。

        通常は直接インスタンス化せず、createファクトリメソッドを使用することを推奨する。

        Args:
            path: タンパク質のファイルパス
            id: タンパク質のID（指定しない場合はファイル名がIDとなる）
        """
        self.__path: Path = Path(path)
        self.__id: str = id if id is not None else self.__path.stem
        self.__domain_events: List[DomainEvent] = []

        # 不変条件の検証
        if not self.__path.exists():
            raise ValueError(f"指定されたパス '{self.__path}' が存在しません")

    @property
    def path(self) -> Path:
        """
        タンパク質のファイルパスを取得する。

        Returns:
            ファイルパス
        """
        return self.__path

    @property
    def id(self) -> str:
        """
        タンパク質のIDを取得する。

        Returns:
            タンパク質のID
        """
        return self.__id

    def __eq__(self, other: object) -> bool:
        """
        等価性比較を行う。エンティティの場合はIDのみで比較する。

        Args:
            other: 比較対象のオブジェクト

        Returns:
            等価であればTrue、そうでなければFalse
        """
        if not isinstance(other, Protein):
            return False
        return self.id == other.id

    def __hash__(self) -> int:
        """
        ハッシュ値を計算する。エンティティの場合はIDのみを使用する。

        Returns:
            ハッシュ値
        """
        return hash(self.id)

    def _register_domain_event(self, event: DomainEvent) -> None:
        """
        ドメインイベントを登録する（内部使用のみ）。

        Args:
            event: 登録するドメインイベント
        """
        self.__domain_events.append(event)

    def release_domain_events(self) -> List[DomainEvent]:
        """
        登録されたドメインイベントを解放し、内部リストをクリアする。

        Returns:
            登録されていたドメインイベントのリスト
        """
        events = self.__domain_events.copy()
        self.__domain_events.clear()
        return events

    def with_new_path(self, new_path: str | Path) -> "Protein":
        """
        新しいパスを持つ新しいProteinインスタンスを作成する。

        Args:
            new_path: 新しいパス

        Returns:
            新しいProteinインスタンス
        """
        old_path = self.path
        new_protein = Protein(path=new_path, id=self.id)

        # パスが変更されたことを表すイベントを登録
        new_protein._register_domain_event(
            ProteinPathUpdated(protein_id=self.id, old_path=str(old_path), new_path=str(new_protein.path))
        )

        return new_protein

    @classmethod
    def create(cls, path: str | Path, id: Optional[str] = None) -> "Protein":
        """
        Proteinオブジェクトを作成するファクトリメソッド。

        Args:
            path: タンパク質のファイルパス
            id: タンパク質のID（指定しない場合はUUIDベースで自動生成）

        Returns:
            作成されたProteinオブジェクト
        """
        # IDが指定されていない場合はUUIDベースで生成
        actual_id = id
        if actual_id is None:
            # UUIDベースのID
            actual_id = f"protein_{uuid.uuid4()}"

        protein = cls(path=path, id=actual_id)

        # タンパク質が登録されたことを表すイベントを登録
        protein._register_domain_event(ProteinRegistered(protein_id=protein.id, path=str(protein.path)))

        return protein

    # TODO: [P2] [DDD] ドメインロジックを追加する
    def get_binding_sites(self) -> List[Dict[str, Any]]:
        """
        タンパク質の結合部位を取得する。

        Returns:
            結合部位のリスト
        """
        # 現時点では空のリストを返す
        # 将来的にはタンパク質の構造解析を行い、結合部位を特定する
        return []

    def get_properties(self) -> Dict[str, Any]:
        """
        タンパク質のプロパティを取得する。

        Returns:
            プロパティの辞書
        """
        # 基本的なプロパティを返す
        return {
            "id": self.id,
            "path": str(self.path),
            "file_format": self.path.suffix.lstrip("."),
            "file_size": self.path.stat().st_size if self.path.exists() else 0,
        }

    def __str__(self) -> str:
        """
        タンパク質の文字列表現を返す。

        Returns:
            タンパク質の文字列表現
        """
        return f"Protein: ID={self.id}, パス={self.path}"

    def __repr__(self) -> str:
        """
        タンパク質の再現可能な文字列表現を返す。

        Returns:
            タンパク質の再現可能な文字列表現
        """
        return f"Protein(path='{self.path}', id='{self.id}')"
