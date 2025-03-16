from dataclasses import dataclass
from typing import Any, Union, Dict, List, Optional, Callable, cast
from enum import Enum, auto


class PropertyType(Enum):
    """分子特性の種類を表す列挙型"""
    MOLECULAR_WEIGHT = auto()
    LOGP = auto()
    NUMBER_OF_ATOMS = auto()
    NUMBER_OF_BONDS = auto()
    NUMBER_OF_ROTATABLE_BONDS = auto()
    NUMBER_OF_HYDROGEN_BOND_DONORS = auto()
    NUMBER_OF_HYDROGEN_BOND_ACCEPTORS = auto()
    POLAR_SURFACE_AREA = auto()
    CHARGE = auto()
    CUSTOM = auto()
    
    @classmethod
    def from_string(cls, name: str) -> 'PropertyType':
        """文字列から特性の種類を取得"""
        name_map = {
            'molecular_weight': cls.MOLECULAR_WEIGHT,
            'logp': cls.LOGP,
            'n_atoms': cls.NUMBER_OF_ATOMS,
            'n_bonds': cls.NUMBER_OF_BONDS,
            'n_rotatable_bonds': cls.NUMBER_OF_ROTATABLE_BONDS,
            'n_hb_donors': cls.NUMBER_OF_HYDROGEN_BOND_DONORS,
            'n_hb_acceptors': cls.NUMBER_OF_HYDROGEN_BOND_ACCEPTORS,
            'polar_surface_area': cls.POLAR_SURFACE_AREA,
            'charge': cls.CHARGE,
        }
        return name_map.get(name.lower(), cls.CUSTOM)


@dataclass(frozen=True)
class MoleculeProperty:
    """分子の物理化学的性質を表す値オブジェクト"""
    name: str
    value: Any
    property_type: PropertyType = PropertyType.CUSTOM
    unit: Optional[str] = None
    description: Optional[str] = None
    
    def __post_init__(self) -> None:
        """初期化後の処理"""
        # もしproperty_typeがCUSTOMで、nameから推測できる場合は設定
        if self.property_type == PropertyType.CUSTOM:
            object.__setattr__(self, 'property_type', PropertyType.from_string(self.name))
    
    def is_valid(self) -> bool:
        """値の妥当性を検証"""
        # 特性タイプごとの値の型チェック
        type_checks = {
            PropertyType.MOLECULAR_WEIGHT: lambda v: isinstance(v, (int, float)) and v > 0,
            PropertyType.LOGP: lambda v: isinstance(v, (int, float)),
            PropertyType.NUMBER_OF_ATOMS: lambda v: isinstance(v, int) and v >= 0,
            PropertyType.NUMBER_OF_BONDS: lambda v: isinstance(v, int) and v >= 0,
            PropertyType.NUMBER_OF_ROTATABLE_BONDS: lambda v: isinstance(v, int) and v >= 0,
            PropertyType.NUMBER_OF_HYDROGEN_BOND_DONORS: lambda v: isinstance(v, int) and v >= 0,
            PropertyType.NUMBER_OF_HYDROGEN_BOND_ACCEPTORS: lambda v: isinstance(v, int) and v >= 0,
            PropertyType.POLAR_SURFACE_AREA: lambda v: isinstance(v, (int, float)) and v >= 0,
            PropertyType.CHARGE: lambda v: isinstance(v, (int, float)),
            PropertyType.CUSTOM: lambda v: True,  # カスタム特性は常に有効とみなす
        }
        
        check_func: Callable[[Any], bool] = cast(Callable[[Any], bool], type_checks.get(self.property_type, lambda v: True))
        return check_func(self.value)
    
    def get_formatted_value(self) -> str:
        """単位付きの値を文字列形式で取得"""
        if self.unit:
            return f"{self.value} {self.unit}"
        return str(self.value)


@dataclass(frozen=True)
class MoleculeProperties:
    """分子の複数の物理化学的特性をまとめた値オブジェクト"""
    properties: List[MoleculeProperty]
    
    def get_property(self, property_type: PropertyType) -> Optional[MoleculeProperty]:
        """指定した種類の特性を取得"""
        for prop in self.properties:
            if prop.property_type == property_type:
                return prop
        return None
    
    def get_property_by_name(self, name: str) -> Optional[MoleculeProperty]:
        """指定した名前の特性を取得"""
        for prop in self.properties:
            if prop.name.lower() == name.lower():
                return prop
        return None
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        return {prop.name: prop.value for prop in self.properties}
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'MoleculeProperties':
        """辞書からMoleculePropertiesを作成"""
        props = [
            MoleculeProperty(name=key, value=value)
            for key, value in data.items()
        ]
        return cls(properties=props)