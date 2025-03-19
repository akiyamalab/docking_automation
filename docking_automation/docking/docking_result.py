from pathlib import Path
from typing import Dict, Any, Optional


# TODO: [DDD] 値オブジェクトとしての実装を強化する
# - 不変性（immutability）を保証する
# - 属性値に基づく等価性（equality）を実装する
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
        id: Optional[str] = None
    ):
        """
        DockingResultオブジェクトを初期化する。
        
        TODO: [DDD] 初期化処理を強化する
        - 一意の識別子（ID）を生成・管理する仕組みを導入する
        - 不変条件のバリデーションを追加する
        - 参照整合性の確認を行う（protein_idとcompound_set_idの存在確認）
        - メタデータのスキーマ検証を追加する
        
        Args:
            result_path: 結果 SDF ファイルのパス
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            docking_score: ドッキングスコア
            metadata: メタデータ
            id: 結果のID（指定しない場合は自動生成）
        """
        self.result_path = result_path
        self.protein_id = protein_id
        self.compound_set_id = compound_set_id
        self.compound_index = compound_index
        self.docking_score = docking_score
        self.metadata = metadata or {}
        self.id = id or f"{protein_id}_{compound_set_id}_{compound_index}"  # TODO: [DDD] より堅牢なID生成方法を実装する
    
    # TODO: [DDD] エンティティの等価性比較を実装する
    def __eq__(self, other: object) -> bool:
        """
        等価性比較を行う。値オブジェクトの場合は全ての属性値で比較する。
        
        Args:
            other: 比較対象のオブジェクト
            
        Returns:
            等価であればTrue、そうでなければFalse
        """
        raise NotImplementedError("値オブジェクトの等価性比較を実装する必要があります")
    
    # TODO: [DDD] 値オブジェクトのハッシュ値計算を実装する
    def __hash__(self) -> int:
        """
        ハッシュ値を計算する。値オブジェクトの場合は全ての属性値を使用する。
        
        Returns:
            ハッシュ値
        """
        raise NotImplementedError("値オブジェクトのハッシュ値計算を実装する必要があります")
    
    # TODO: [DDD] ドメインロジックを追加する
    def get_pose(self) -> Path:
        """
        ドッキングポーズのファイルパスを取得する。
        
        Returns:
            ポーズファイルのパス
        """
        raise NotImplementedError("ポーズ取得機能を実装する必要があります")
    
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
    
    # TODO: [DDD] 値オブジェクトの不変性を保持するメソッドを実装する
    def with_metadata(self, key: str, value: Any) -> 'DockingResult':
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
            metadata=new_metadata
        )
    
    # TODO: [DDD] 値オブジェクト生成のためのファクトリメソッドを追加する
    @classmethod
    def create(
        cls,
        result_path: Path,
        protein_id: str,
        compound_set_id: str,
        compound_index: int,
        docking_score: float,
        metadata: Optional[Dict[str, Any]] = None
    ) -> 'DockingResult':
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
            
        Returns:
            作成されたDockingResultオブジェクト
        """
        raise NotImplementedError("値オブジェクト生成のファクトリメソッドを実装する必要があります")
    
    # TODO: [DDD] 値オブジェクトの変更を通知する機能を追加する
    def with_domain_event(self, event: Any) -> 'DockingResult':  # 'DomainEvent'型は将来実装予定
        """
        ドメインイベントを含む新しいDockingResultインスタンスを作成する。
        値オブジェクトは不変なので、イベントを登録した新しいインスタンスを返す。
        
        Args:
            event: 登録するドメインイベント
            
        Returns:
            ドメインイベントを含む新しいDockingResultインスタンス
        """
        raise NotImplementedError("値オブジェクトの変更通知機能を実装する必要があります")