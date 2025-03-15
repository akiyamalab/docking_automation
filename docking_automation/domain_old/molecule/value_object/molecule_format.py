from enum import Enum, auto
from dataclasses import dataclass
from typing import Optional


class FormatType(Enum):
    """分子形式の種類を表す列挙型"""
    SDF = auto()
    PDB = auto()
    PDBQT = auto()
    MOL2 = auto()
    SMILES = auto()
    UNKNOWN = auto()
    
    @classmethod
    def from_extension(cls, file_extension: str) -> 'FormatType':
        """ファイル拡張子から形式を推測"""
        extension_map = {
            'sdf': cls.SDF,
            'mol': cls.SDF,
            'pdb': cls.PDB,
            'pdbqt': cls.PDBQT,
            'mol2': cls.MOL2,
            'smi': cls.SMILES,
            'smiles': cls.SMILES,
        }
        return extension_map.get(file_extension.lower().lstrip('.'), cls.UNKNOWN)
    
    @classmethod
    def from_path(cls, file_path: str) -> 'FormatType':
        """ファイルパスから形式を推測"""
        if not file_path:
            return cls.UNKNOWN
        
        extension = file_path.split('.')[-1] if '.' in file_path else ''
        return cls.from_extension(extension)
    
    def get_extension(self) -> str:
        """対応するファイル拡張子を取得"""
        if self == FormatType.SDF:
            return 'sdf'
        elif self == FormatType.PDB:
            return 'pdb'
        elif self == FormatType.PDBQT:
            return 'pdbqt'
        elif self == FormatType.MOL2:
            return 'mol2'
        elif self == FormatType.SMILES:
            return 'smi'
        else:
            return 'txt'


@dataclass(frozen=True)
class MoleculeFormat:
    """分子形式を表す値オブジェクト"""
    type: FormatType
    version: Optional[str] = None
    
    def validate(self) -> bool:
        """形式の妥当性を検証"""
        return self.type != FormatType.UNKNOWN
    
    def is_3d_format(self) -> bool:
        """3D座標を含む形式かどうかを判定"""
        return self.type in [FormatType.SDF, FormatType.PDB, FormatType.PDBQT, FormatType.MOL2]
    
    def is_compatible_with(self, other: 'MoleculeFormat') -> bool:
        """他の形式と互換性があるかどうかを判定"""
        # 同じ形式は互換性あり
        if self.type == other.type:
            return True
        
        # 3D形式同士は基本的に変換可能
        if self.is_3d_format() and other.is_3d_format():
            return True
        
        # SMILESから3D形式への変換は可能だが、情報が失われる可能性がある
        if self.type == FormatType.SMILES and other.is_3d_format():
            return False
        
        return False