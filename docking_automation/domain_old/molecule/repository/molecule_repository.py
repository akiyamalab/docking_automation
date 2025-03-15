from abc import ABC, abstractmethod
from typing import List, Optional, Generic, TypeVar, Protocol, Dict, Any

from docking_automation.domain.molecule.entity.molecule import Molecule
from docking_automation.domain.molecule.value_object.molecule_format import FormatType


# 型変数の定義: T は Molecule を継承した任意の型
T = TypeVar('T', bound=Molecule)


class MoleculeRepository(Generic[T], ABC):
    """分子リポジトリの抽象基底クラス
    
    このリポジトリは、分子エンティティの保存と取得を担当します。
    """
    
    @abstractmethod
    def save(self, molecule: T) -> T:
        """分子を保存する
        
        Args:
            molecule: 保存する分子
            
        Returns:
            保存された分子（IDが更新される可能性あり）
            
        Raises:
            IOError: 保存処理中にエラーが発生した場合
        """
        pass
    
    @abstractmethod
    def find_by_id(self, id: str) -> Optional[T]:
        """IDで分子を検索する
        
        Args:
            id: 検索する分子のID
            
        Returns:
            見つかった分子、見つからない場合はNone
            
        Raises:
            IOError: 検索処理中にエラーが発生した場合
        """
        pass
    
    @abstractmethod
    def find_all(self) -> List[T]:
        """すべての分子を取得する
        
        Returns:
            すべての分子のリスト
            
        Raises:
            IOError: 取得処理中にエラーが発生した場合
        """
        pass
    
    @abstractmethod
    def find_by_property(self, property_name: str, property_value: Any) -> List[T]:
        """特定のプロパティ値を持つ分子を検索する
        
        Args:
            property_name: 検索するプロパティ名
            property_value: 検索するプロパティ値
            
        Returns:
            条件に一致する分子のリスト
            
        Raises:
            IOError: 検索処理中にエラーが発生した場合
        """
        pass
    
    @abstractmethod
    def delete(self, id: str) -> bool:
        """指定されたIDの分子を削除する
        
        Args:
            id: 削除する分子のID
            
        Returns:
            削除が成功したかどうか
            
        Raises:
            IOError: 削除処理中にエラーが発生した場合
        """
        pass
    
    @abstractmethod
    def load_from_file(self, file_path: str, format_type: Optional[FormatType] = None) -> T:
        """ファイルから分子をロードする
        
        Args:
            file_path: ロードするファイルのパス
            format_type: ファイル形式（Noneの場合は拡張子から推測）
            
        Returns:
            ロードされた分子
            
        Raises:
            IOError: ファイル読み込みに失敗した場合
            ValueError: ファイル形式が無効な場合
        """
        pass
    
    @abstractmethod
    def save_to_file(self, molecule: T, file_path: str, format_type: Optional[FormatType] = None) -> bool:
        """分子をファイルに保存する
        
        Args:
            molecule: 保存する分子
            file_path: 保存先ファイルのパス
            format_type: ファイル形式（Noneの場合は拡張子から推測）
            
        Returns:
            保存が成功したかどうか
            
        Raises:
            IOError: ファイル書き込みに失敗した場合
            ValueError: ファイル形式が無効な場合
        """
        pass
    
    @abstractmethod
    def query(self, query_params: Dict[str, Any]) -> List[T]:
        """高度な検索条件で分子を検索する
        
        Args:
            query_params: 検索パラメータの辞書
            
        Returns:
            条件に一致する分子のリスト
            
        Raises:
            IOError: 検索処理中にエラーが発生した場合
            ValueError: 無効な検索パラメータが指定された場合
        """
        pass