import hashlib
import os
import shutil
import tempfile
import uuid
from functools import cached_property
from pathlib import Path
from typing import Any, Dict, Iterator, List, Optional, Tuple

from docking_automation.domain.domain_event import DomainEvent
from docking_automation.infrastructure.utilities.file_utils import (
    calculate_file_content_hash,
    count_compounds_in_sdf,
    read_compounds_from_sdf,
    safe_open,
)
from docking_automation.molecule.compound_set_events import (
    CompoundSetPathUpdated,
    CompoundSetRegistered,
)

# TODO: [P1] [DDD] インフラストラクチャの依存性（RDKit）をドメインモデルから分離する
# TODO: [P1] [DDD] リポジトリパターンを導入し、CompoundSetエンティティの永続化を担当するリポジトリを実装する
# TODO: [P2] [DDD] 集約の境界を明確にする（CompoundSetは複数のCompoundエンティティの集約ルートとして機能すべき）


class CompoundSet:
    """
    化合物セットを表すクラス。
    現時点では .sdf 形式をサポートしている。
    1つのファイルに複数の化合物が含まれる。

    IDを持つ不変エンティティとして実装されており、初期化後に属性を変更することはできない。
    状態変化が必要な場合は、新しいインスタンスを作成する。
    """

    def __init__(self, path: str | Path, id: Optional[str] = None):
        """
        CompoundSetオブジェクトを初期化する。

        通常は直接インスタンス化せず、createファクトリメソッドを使用することを推奨する。

        Args:
            path: 化合物セットのファイルパス
            id: 化合物セットのID（指定しない場合はファイル名がIDとなる）
        """
        self.__path: Path = Path(path)
        self.__id: str = id if id is not None else self.__path.stem
        self.__domain_events: List[DomainEvent] = []
        self.__index_range: Optional[Tuple[int, int]] = None  # インデックス範囲（開始、終了）
        self.__compound_count: Optional[int] = None  # 化合物数のキャッシュ
        self.__compound_hash_cache: Dict[int, str] = {}  # インデックス -> 化合物ハッシュ値のキャッシュ

        # 一時ファイル関連の属性（デフォルトではNone）
        # これらの属性は分割時に設定される
        self._original_indices: Optional[List[int]] = None  # 元のインデックスのリスト
        self._original_compound_set_id: Optional[str] = None  # 元の化合物セットのID
        self._temp_dir: Optional[str] = None  # 一時ディレクトリのパス

        # 不変条件の検証
        if not self.__path.exists():
            raise ValueError(f"指定されたパス '{self.__path}' が存在しません")

        # 初期化時にすべての化合物のハッシュ値を計算
        self._calculate_all_compound_hashes()

    @property
    def path(self) -> Path:
        """
        化合物セットのファイルパスを取得する。

        Returns:
            ファイルパス
        """
        return self.__path

    @property
    def id(self) -> str:
        """
        化合物セットのIDを取得する。

        Returns:
            化合物セットのID
        """
        return self.__id

    def get_compound_count(self) -> int:
        """
        化合物セットに含まれる化合物の数を取得する。

        パスが更新されたときにのみ実際のカウントを行い、それ以外の場合はキャッシュされた値を返す。

        TODO: [P1] [DDD] インフラストラクチャの依存性を分離する
        - ファイル読み込みロジックをリポジトリに委譲する

        Returns:
            化合物の数
        """
        # キャッシュされた値がある場合はそれを返す
        if self.__compound_count is not None:
            return self.__compound_count

        # キャッシュがない場合は実際にカウントする
        # ファイル形式を判断
        file_format = self.path.suffix.lower()

        # gzipで圧縮されている場合
        if file_format == ".gz":
            # 実際のファイル形式を取得
            actual_format = self.path.stem.split(".")[-1].lower()

            if actual_format == "sdf":
                self.__compound_count = count_compounds_in_sdf(self.path)
            else:
                raise ValueError(f"サポートされていないファイル形式です: {actual_format}")
        # 非圧縮ファイルの場合
        elif file_format == ".sdf":
            self.__compound_count = count_compounds_in_sdf(self.path)
        else:
            raise ValueError(f"サポートされていないファイル形式です: {file_format}")

        return self.__compound_count

    def __eq__(self, other: object) -> bool:
        """
        等価性比較を行う。エンティティの場合はIDのみで比較する。

        Args:
            other: 比較対象のオブジェクト

        Returns:
            等価であればTrue、そうでなければFalse
        """
        if not isinstance(other, CompoundSet):
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

    def with_new_path(self, new_path: str | Path) -> "CompoundSet":
        """
        新しいパスを持つ新しいCompoundSetインスタンスを作成する。

        Args:
            new_path: 新しいパス

        Returns:
            新しいCompoundSetインスタンス
        """
        old_path = self.path
        new_compound_set = CompoundSet(path=new_path, id=self.id)

        # 新しいパスでは化合物数を再計算する必要があるため、キャッシュを無効化
        # 次回get_compound_countが呼ばれたときに再計算される
        new_compound_set.__compound_count = None

        # パスが変更されたことを表すイベントを登録
        new_compound_set._register_domain_event(
            CompoundSetPathUpdated(compound_set_id=self.id, old_path=str(old_path), new_path=str(new_compound_set.path))
        )

        return new_compound_set

    def split_by_chunks(self, chunk_size: int, compress: bool = False) -> List["CompoundSet"]:
        """
        化合物セットを指定されたチャンクサイズで分割する。

        大きな化合物セットを小さなセットに分割して、並列処理などを効率化するために使用する。
        分割された各化合物セットは元のファイルの一部を含む新しい一時ファイルを参照する。
        これにより、複数のタスクが同時に同じファイルにアクセスする際の非効率性を回避する。
        分割されたファイルは一時ファイルとして保存され、計算終了後に自動的に削除される。

        Args:
            chunk_size: 各チャンクに含める化合物の最大数
            compress: 分割後のファイルを圧縮するかどうか（デフォルトはFalse）

        Returns:
            分割された化合物セットのリスト

        Raises:
            ValueError: チャンクサイズが1未満の場合
        """
        if chunk_size < 1:
            raise ValueError("チャンクサイズは1以上である必要があります")

        # 分割されたCompoundSetのリスト
        compound_sets = []
        # 開始インデックス
        start_index = 0
        # 化合物の総数
        total_compounds = self.get_compound_count()

        # 各チャンクごとに独自の一時ディレクトリを作成
        # これにより、各CompoundSetオブジェクトが独立した一時ディレクトリを持つようになる
        main_temp_dir = tempfile.mkdtemp()

        # ファイル形式を判断

        try:
            # SDFファイルから化合物データを読み込む
            # TODO: このやり方だと一旦ファイルを全て読み込み保存する必要があり、メモリを圧迫する可能性がある
            #       代わりに、ファイルをストリームとして読み込む方法を検討すべき
            compound_data = list(read_compounds_from_sdf(self.path))

            # チャンクごとに一時ファイルを作成（各チャンクは独自の一時ディレクトリを持つ）
            for i in range(0, total_compounds, chunk_size):
                end_index = min(i + chunk_size, total_compounds)
                chunk_compounds = compound_data[i:end_index]

                # チャンクごとに独自の一時ディレクトリを作成
                chunk_temp_dir = os.path.join(main_temp_dir, f"chunk_{i}_{end_index}")
                os.makedirs(chunk_temp_dir, exist_ok=True)

                # 一時ファイルのパスを生成
                temp_file_path = os.path.join(chunk_temp_dir, f"{self.id}_chunk_{i}_{end_index}.sdf")
                # 圧縮オプションに基づいてファイル形式を決定
                if compress:
                    temp_file_path += ".gz"

                # 一時ファイルに書き込み（file_utilsのsafe_open関数を使用）
                with safe_open(temp_file_path, "wt") as temp_f:
                    for _, compound_lines in chunk_compounds:
                        for line in compound_lines:
                            temp_f.write(line)  # type: ignore

                # 新しいCompoundSetを作成し、元のインデックス情報を保持
                new_compound_set = CompoundSet(path=temp_file_path, id=f"{self.id}_chunk_{i}_{end_index}")

                # 元のインデックス情報をメタデータとして保存
                original_indices = [idx for idx, _ in chunk_compounds]
                new_compound_set._original_indices = original_indices
                new_compound_set._original_compound_set_id = self.id
                new_compound_set._temp_dir = str(chunk_temp_dir)  # チャンク固有の一時ディレクトリへの参照を保持

                compound_sets.append(new_compound_set)

            return compound_sets
        except Exception as e:
            # エラーが発生した場合はメイン一時ディレクトリを削除
            shutil.rmtree(main_temp_dir, ignore_errors=True)
            raise ValueError(f"化合物セットの分割中にエラーが発生しました: {e}")

    def with_index_range(self, start_index: int, end_index: int) -> "CompoundSet":
        """
        指定されたインデックス範囲の化合物のみを含む新しいCompoundSetを作成する。

        元のCompoundSetと同じファイルを参照するが、処理対象のインデックス範囲が制限される。
        不変性を維持するため、新しいインスタンスを返す。

        Args:
            start_index: 開始インデックス（含む）
            end_index: 終了インデックス（含まない）

        Returns:
            インデックス範囲が制限された新しいCompoundSet

        Raises:
            ValueError: インデックスが範囲外の場合
        """
        # インデックス範囲の検証
        total_compounds = self.get_compound_count()
        if start_index < 0 or start_index >= total_compounds:
            raise ValueError(f"開始インデックス {start_index} は範囲外です（0-{total_compounds-1}）")
        if end_index <= start_index or end_index > total_compounds:
            raise ValueError(f"終了インデックス {end_index} は無効です（{start_index+1}-{total_compounds}）")

        # 新しいCompoundSetを作成（同じファイルを参照）
        new_compound_set = CompoundSet(path=self.path, id=f"{self.id}_range_{start_index}_{end_index}")

        # インデックス範囲を保存するための属性を追加
        # プライベート属性として追加し、外部からは直接アクセスできないようにする
        new_compound_set.__index_range = (start_index, end_index)

        # ドメインイベントを登録（オプション）
        # 新しいCompoundSetが作成されたことを表すイベントを登録
        # 現時点では特定のイベントクラスがないため、汎用的なイベントを使用

        return new_compound_set

    # TODO: [P2] [DDD] 集約ルートとしての操作を追加する
    def get_compound(self, index: int) -> Dict[str, Any]:
        """
        指定されたインデックスの化合物を取得する。

        インデックス範囲が設定されている場合は、その範囲内のインデックスのみが有効。
        分割されたCompoundSetの場合、元のCompoundSetにおけるインデックスも返される。

        Args:
            index: 取得する化合物のインデックス

        Returns:
            化合物の情報を含む辞書

        Raises:
            IndexError: インデックスが範囲外の場合
        """
        # インデックス範囲が設定されている場合は、その範囲内かどうかをチェック
        if self.__index_range is not None:
            start_index, end_index = self.__index_range
            if index < start_index or index >= end_index:
                raise IndexError(f"インデックス {index} は設定された範囲外です（{start_index}-{end_index-1}）")

        # 全体の範囲内かどうかをチェック
        total_compounds = self.get_compound_count()
        if index < 0 or index >= total_compounds:
            raise IndexError(f"インデックス {index} は範囲外です（0-{total_compounds-1}）")

        # 基本情報
        compound_info = {"index": index, "compound_set_id": self.id}

        # 元のインデックス情報がある場合は追加
        if (
            hasattr(self, "_original_indices")
            and self._original_indices is not None
            and index < len(self._original_indices)
        ):
            compound_info["original_index"] = self._original_indices[index]
            compound_info["original_compound_set_id"] = self._original_compound_set_id

        return compound_info

    # TODO: [P2] [DDD] 集約ルートとしての操作を追加する
    def get_all_compounds(self) -> List[Dict[str, Any]]:
        """
        すべての化合物を取得する。

        インデックス範囲が設定されている場合は、その範囲内の化合物のみを返す。

        Returns:
            化合物のリスト
        """
        # インデックス範囲が設定されている場合は、その範囲内の化合物のみを返す
        if self.__index_range is not None:
            start_index, end_index = self.__index_range
            return [self.get_compound(i) for i in range(start_index, end_index)]

        # インデックス範囲が設定されていない場合は、すべての化合物を返す
        return [self.get_compound(i) for i in range(self.get_compound_count())]

    @classmethod
    def create(cls, path: str | Path, id: Optional[str] = None) -> "CompoundSet":
        """
        CompoundSetオブジェクトを作成するファクトリメソッド。

        Args:
            path: 化合物セットのファイルパス
            id: 化合物セットのID（指定しない場合はUUIDベースで自動生成）

        Returns:
            作成されたCompoundSetオブジェクト
        """
        # IDが指定されていない場合はUUIDベースで生成
        actual_id = id
        if actual_id is None:
            path_obj = Path(path)
            if path_obj.exists():
                # ファイル名をベースにしたID
                actual_id = path_obj.stem
            else:
                # UUIDベースのID
                actual_id = f"compound_set_{uuid.uuid4()}"

        compound_set = cls(path=path, id=actual_id)

        # 化合物セットが登録されたことを表すイベントを登録
        compound_set._register_domain_event(
            CompoundSetRegistered(compound_set_id=compound_set.id, path=str(compound_set.path))
        )

        return compound_set

    @classmethod
    def create_empty(
        cls, id: Optional[str] = None, temp_dir: Optional[str | Path] = None, compress: bool = False
    ) -> "CompoundSet":
        """
        空の化合物セットを作成するファクトリメソッド。

        Args:
            id: 化合物セットのID（指定しない場合は自動生成）
            temp_dir: 一時ファイルを作成するディレクトリ
            compress: ファイルを圧縮するかどうか（デフォルトはFalse）

        Returns:
            作成された空の化合物セット
        """

        # 一時ファイルを作成
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()

        temp_dir_path = Path(temp_dir)
        os.makedirs(temp_dir_path, exist_ok=True)

        # UUIDベースのIDを生成
        actual_id = id or f"compound_set_{uuid.uuid4()}"

        # 空のSDFファイルを作成（圧縮オプションに基づいてファイル形式を決定）
        temp_file = temp_dir_path / f"{actual_id}.sdf"
        if compress:
            temp_file = temp_file.with_suffix(".sdf.gz")

        with safe_open(temp_file, "w") as f:
            # 空のSDFファイル
            pass

        # CompoundSetオブジェクトを作成
        compound_set = cls.create(path=temp_file, id=actual_id)
        compound_set._temp_dir = str(temp_dir) if temp_dir is not None else None  # 一時ディレクトリへの参照を保持
        return compound_set

    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """
        化合物セットのイテレータを取得する。

        インデックス範囲が設定されている場合は、その範囲内の化合物のみをイテレートする。
        分割されたCompoundSetの場合、元のCompoundSetにおけるインデックス情報も含まれる。

        Returns:
            化合物のイテレータ
        """
        # インデックス範囲が設定されている場合は、その範囲内の化合物のみをイテレートする
        if self.__index_range is not None:
            start_index, end_index = self.__index_range
            for i in range(start_index, end_index):
                yield self.get_compound(i)
        else:
            # インデックス範囲が設定されていない場合は、すべての化合物をイテレートする
            for i in range(self.get_compound_count()):
                yield self.get_compound(i)

    def __len__(self) -> int:
        """
        化合物セットの長さ（化合物の数）を取得する。

        インデックス範囲が設定されている場合は、その範囲内の化合物の数を返す。

        Returns:
            化合物の数
        """
        # インデックス範囲が設定されている場合は、その範囲内の化合物の数を返す
        if self.__index_range is not None:
            start_index, end_index = self.__index_range
            return end_index - start_index

        # インデックス範囲が設定されていない場合は、すべての化合物の数を返す
        return self.get_compound_count()

    def get_properties(self) -> Dict[str, Any]:
        """
        化合物セットのプロパティを取得する。

        インデックス範囲が設定されている場合は、その情報も含める。
        分割されたCompoundSetの場合、元のCompoundSetに関する情報も含める。

        Returns:
            プロパティの辞書
        """
        # 基本的なプロパティを返す
        properties = {
            "id": self.id,
            "path": str(self.path),
            "file_format": self.path.suffix.lstrip("."),
            "file_size": self.path.stat().st_size if self.path.exists() else 0,
            "compound_count": len(
                self
            ),  # __len__メソッドを使用して、インデックス範囲が設定されている場合はその範囲内の化合物の数を返す
            "content_hash": self.content_hash,  # ファイル内容のハッシュ値
        }

        # インデックス範囲が設定されている場合は、その情報も含める
        if self.__index_range is not None:
            start_index, end_index = self.__index_range
            properties["index_range"] = {
                "start": start_index,
                "end": end_index,
                "total_compounds": self.get_compound_count(),  # 全体の化合物の数
            }

        # 元のインデックス情報がある場合は追加
        if hasattr(self, "_original_indices"):
            properties["original_indices"] = self._original_indices
            properties["original_compound_set_id"] = self._original_compound_set_id

        return properties

    def _calculate_all_compound_hashes(self) -> None:
        """すべての化合物のハッシュ値を計算してキャッシュする"""
        for i, (idx, compound_lines) in enumerate(read_compounds_from_sdf(self.path)):
            # 化合物のテキスト表現からハッシュ値を計算
            content = "".join(compound_lines)
            hash_value = hashlib.sha256(content.encode("utf-8")).hexdigest()
            self.__compound_hash_cache[i] = hash_value

    def get_compound_hash(self, index: int) -> str:
        """
        指定されたインデックスの化合物のハッシュ値を取得する

        Args:
            index: 取得する化合物のインデックス

        Returns:
            化合物のハッシュ値（SHA-256ハッシュ値の16進数文字列）

        Raises:
            IndexError: インデックスが範囲外の場合
        """
        if index not in self.__compound_hash_cache:
            raise IndexError(f"インデックス {index} は範囲外です")
        return self.__compound_hash_cache[index]

    @cached_property
    def content_hash(self) -> str:
        """
        ファイルの内容に基づいたハッシュ値を取得する。

        初回アクセス時に計算し、その後はキャッシュした値を返す。

        Returns:
            ファイル内容のSHA-256ハッシュ値（16進数文字列）
        """
        return calculate_file_content_hash(self.path)
