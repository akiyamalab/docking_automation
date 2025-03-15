from abc import ABC, abstractmethod
from typing import TypeVar, Generic, Optional, List, Tuple

from docking_automation.molecule.entity.molecule import Molecule
from docking_automation.molecule.value_object.molecule_format import FormatType


# 型変数の定義: MoleculeTypeは任意のMoleculeサブクラス
MoleculeType = TypeVar('MoleculeType', bound=Molecule)


class FormatConverterService(Generic[MoleculeType], ABC):
    """分子形式変換サービスのインターフェース
    
    このサービスは、分子の形式（例：SDF, PDB, PDBQT）間の変換を担当します。
    """
    
    @abstractmethod
    def convert(self, molecule: MoleculeType, target_format: FormatType) -> MoleculeType:
        """分子を指定された形式に変換する
        
        Args:
            molecule: 変換する分子
            target_format: 変換先の形式
            
        Returns:
            変換された分子（元の分子と同じ型）
        
        Raises:
            ValueError: 変換元と変換先の形式が互換性がない場合
            IOError: ファイル操作に失敗した場合
        """
        pass
    
    @abstractmethod
    def convert_file(self, input_path: str, output_path: str, target_format: FormatType) -> bool:
        """ファイルを指定された形式に変換する
        
        Args:
            input_path: 入力ファイルのパス
            output_path: 出力ファイルのパス
            target_format: 変換先の形式
            
        Returns:
            変換が成功したかどうか
            
        Raises:
            ValueError: 入力ファイルが無効な場合
            IOError: ファイル操作に失敗した場合
        """
        pass
    
    @abstractmethod
    def get_supported_conversions(self) -> List[Tuple[FormatType, FormatType]]:
        """このサービスがサポートする変換のリストを取得する
        
        Returns:
            サポートされている変換のリスト（入力形式、出力形式）のタプル
        """
        pass
    
    @abstractmethod
    def can_convert(self, source_format: FormatType, target_format: FormatType) -> bool:
        """指定された変換がサポートされているかどうかを判定する
        
        Args:
            source_format: 変換元の形式
            target_format: 変換先の形式
            
        Returns:
            変換がサポートされているかどうか
        """
        pass