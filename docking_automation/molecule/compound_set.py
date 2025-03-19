import gzip
import uuid
from pathlib import Path
from typing import Optional, List, Dict, Any, Iterator, Callable, Set

from docking_automation.domain.domain_event import DomainEvent
from docking_automation.molecule.compound_set_events import CompoundSetRegistered, CompoundSetPathUpdated

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
    def __init__(self, path: str|Path, id: Optional[str] = None):
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
        
        # 不変条件の検証
        if not self.__path.exists():
            raise ValueError(f"指定されたパス '{self.__path}' が存在しません")
    
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
        
        TODO: [P1] [DDD] インフラストラクチャの依存性を分離する
        - ファイル読み込みロジックをリポジトリに委譲する
        
        Returns:
            化合物の数
        """
        # ファイル形式を判断
        file_format = self.path.suffix.lower()
        
        # gzipで圧縮されている場合
        if file_format == '.gz':
            # 実際のファイル形式を取得
            actual_format = self.path.stem.split('.')[-1].lower()
            
            if actual_format == 'sdf':
                return self._count_compounds_in_sdf()
            else:
                raise ValueError(f"サポートされていないファイル形式です: {actual_format}")
        # 非圧縮ファイルの場合
        elif file_format == '.sdf':
            return self._count_compounds_in_sdf()
        else:
            raise ValueError(f"サポートされていないファイル形式です: {file_format}")
    
    def _count_compounds_in_sdf(self) -> int:
        """
        SDFファイル内の化合物数をカウントする。
        
        Returns:
            化合物の数
        """
        try:
            if str(self.path).endswith('.gz'):
                with gzip.open(self.path, 'rt') as f:
                    count = 0
                    for line in f:
                        if "$$$$" in line:
                            count += 1
                    return count
            else:
                with open(self.path, 'rt') as f:
                    count = 0
                    for line in f:
                        if "$$$$" in line:
                            count += 1
                    return count
        except Exception as e:
            raise ValueError(f"化合物数のカウント中にエラーが発生しました: {e}")
    
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
    
    def with_new_path(self, new_path: str|Path) -> 'CompoundSet':
        """
        新しいパスを持つ新しいCompoundSetインスタンスを作成する。
        
        Args:
            new_path: 新しいパス
            
        Returns:
            新しいCompoundSetインスタンス
        """
        old_path = self.path
        new_compound_set = CompoundSet(path=new_path, id=self.id)
        
        # パスが変更されたことを表すイベントを登録
        new_compound_set._register_domain_event(CompoundSetPathUpdated(
            compound_set_id=self.id,
            old_path=str(old_path),
            new_path=str(new_compound_set.path)
        ))
        
        return new_compound_set
    
    # TODO: [P2] [DDD] 集約ルートとしての操作を追加する
    def get_compound(self, index: int) -> Dict[str, Any]:
        """
        指定されたインデックスの化合物を取得する。
        
        Args:
            index: 取得する化合物のインデックス
            
        Returns:
            化合物の情報を含む辞書
        """
        # 現時点では空の辞書を返す
        # 将来的には実際の化合物情報を返す
        if index < 0 or index >= self.get_compound_count():
            raise IndexError(f"インデックス {index} は範囲外です（0-{self.get_compound_count()-1}）")
        
        return {
            "index": index,
            "compound_set_id": self.id
        }
    
    # TODO: [P2] [DDD] 集約ルートとしての操作を追加する
    def get_all_compounds(self) -> List[Dict[str, Any]]:
        """
        すべての化合物を取得する。
        
        Returns:
            化合物のリスト
        """
        # 現時点では空のリストを返す
        # 将来的には実際の化合物情報のリストを返す
        return [self.get_compound(i) for i in range(self.get_compound_count())]
    
    @classmethod
    def create(cls, path: str|Path, id: Optional[str] = None) -> 'CompoundSet':
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
        compound_set._register_domain_event(CompoundSetRegistered(
            compound_set_id=compound_set.id,
            path=str(compound_set.path)
        ))
        
        return compound_set
    
    @classmethod
    def create_empty(cls, id: Optional[str] = None, temp_dir: Optional[str|Path] = None) -> 'CompoundSet':
        """
        空の化合物セットを作成するファクトリメソッド。
        
        Args:
            id: 化合物セットのID（指定しない場合は自動生成）
            temp_dir: 一時ファイルを作成するディレクトリ
            
        Returns:
            作成された空の化合物セット
        """
        import tempfile
        import os
        
        # 一時ファイルを作成
        if temp_dir is None:
            temp_dir = tempfile.gettempdir()
        
        temp_dir_path = Path(temp_dir)
        os.makedirs(temp_dir_path, exist_ok=True)
        
        # UUIDベースのIDを生成
        actual_id = id or f"compound_set_{uuid.uuid4()}"
        
        # 空のSDFファイルを作成
        temp_file = temp_dir_path / f"{actual_id}.sdf"
        with open(temp_file, 'w') as f:
            # 空のSDFファイル
            pass
        
        # CompoundSetオブジェクトを作成
        return cls.create(path=temp_file, id=actual_id)
    
    def __iter__(self) -> Iterator[Dict[str, Any]]:
        """
        化合物セットのイテレータを取得する。
        
        Returns:
            化合物のイテレータ
        """
        for i in range(self.get_compound_count()):
            yield self.get_compound(i)
    
    def __len__(self) -> int:
        """
        化合物セットの長さ（化合物の数）を取得する。
        
        Returns:
            化合物の数
        """
        return self.get_compound_count()
    
    def get_properties(self) -> Dict[str, Any]:
        """
        化合物セットのプロパティを取得する。
        
        Returns:
            プロパティの辞書
        """
        # 基本的なプロパティを返す
        return {
            "id": self.id,
            "path": str(self.path),
            "file_format": self.path.suffix.lstrip('.'),
            "file_size": self.path.stat().st_size if self.path.exists() else 0,
            "compound_count": self.get_compound_count()
        }