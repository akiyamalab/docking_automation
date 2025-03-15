from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Any, Dict, Optional, Union


class ParameterType(Enum):
    """ドッキングパラメータの型を表す列挙型"""
    INTEGER = auto()
    FLOAT = auto()
    BOOLEAN = auto()
    STRING = auto()
    
    @classmethod
    def from_value(cls, value: Any) -> 'ParameterType':
        """値の型からパラメータ型を推測"""
        if isinstance(value, bool):
            return cls.BOOLEAN
        elif isinstance(value, int):
            return cls.INTEGER
        elif isinstance(value, float):
            return cls.FLOAT
        else:
            return cls.STRING


@dataclass(frozen=True)
class DockingParameter:
    """ドッキングパラメータを表す値オブジェクト"""
    name: str
    value: Any
    parameter_type: ParameterType
    description: Optional[str] = None
    min_value: Optional[Union[int, float]] = None
    max_value: Optional[Union[int, float]] = None
    
    def __post_init__(self) -> None:
        """初期化後の検証"""
        # 値の型チェック
        
        # 値の型チェック
        if self.parameter_type == ParameterType.INTEGER and not isinstance(self.value, int):
            raise ValueError(f"Parameter {self.name} must be an integer")
        if self.parameter_type == ParameterType.FLOAT and not isinstance(self.value, (int, float)):
            raise ValueError(f"Parameter {self.name} must be a float")
        if self.parameter_type == ParameterType.BOOLEAN and not isinstance(self.value, bool):
            raise ValueError(f"Parameter {self.name} must be a boolean")
        
        # 範囲チェック（数値パラメータの場合）
        if (self.parameter_type in [ParameterType.INTEGER, ParameterType.FLOAT]
                and (self.min_value is not None or self.max_value is not None)):
            
            if self.min_value is not None and self.value < self.min_value:
                raise ValueError(f"Parameter {self.name} must be >= {self.min_value}")
            
            if self.max_value is not None and self.value > self.max_value:
                raise ValueError(f"Parameter {self.name} must be <= {self.max_value}")
    
    def is_valid(self) -> bool:
        """パラメータが有効かどうかを検証"""
        # 値の型チェック
        if self.parameter_type == ParameterType.INTEGER and not isinstance(self.value, int):
            return False
        if self.parameter_type == ParameterType.FLOAT and not isinstance(self.value, (int, float)):
            return False
        if self.parameter_type == ParameterType.BOOLEAN and not isinstance(self.value, bool):
            return False
        
        # 範囲チェック（数値パラメータの場合）
        if (self.parameter_type in [ParameterType.INTEGER, ParameterType.FLOAT]
                and (self.min_value is not None or self.max_value is not None)):
            
            if self.min_value is not None and self.value < self.min_value:
                return False
            
            if self.max_value is not None and self.value > self.max_value:
                return False
                
        return True
    
    def to_str(self) -> str:
        """パラメータを文字列形式に変換"""
        if self.parameter_type == ParameterType.BOOLEAN:
            return str(self.value).lower()
        return str(self.value)
    
    def to_dict(self) -> Dict[str, Any]:
        """パラメータを辞書形式に変換"""
        result = {
            'name': self.name,
            'value': self.value,
            'type': self.parameter_type.name,
        }
        
        if self.description:
            result['description'] = self.description
        
        if self.min_value is not None:
            result['min_value'] = self.min_value
        
        if self.max_value is not None:
            result['max_value'] = self.max_value
        
        return result


@dataclass
class DockingParameters:
    """ドッキングパラメータのコレクションを表す値オブジェクト"""
    parameters: Dict[str, DockingParameter] = field(default_factory=dict)
    
    def add(self, parameter: DockingParameter) -> None:
        """パラメータを追加"""
        self.parameters[parameter.name] = parameter
    
    def get(self, name: str) -> Optional[DockingParameter]:
        """名前でパラメータを取得"""
        return self.parameters.get(name)
    
    def get_value(self, name: str, default: Any = None) -> Any:
        """名前でパラメータの値を取得"""
        param = self.get(name)
        if param is None:
            return default
        return param.value
    
    def remove(self, name: str) -> None:
        """名前でパラメータを削除"""
        if name in self.parameters:
            del self.parameters[name]
    
    def has(self, name: str) -> bool:
        """指定された名前のパラメータが存在するかどうかを確認"""
        return name in self.parameters
    
    def validate(self) -> bool:
        """すべてのパラメータが有効かどうかを検証"""
        return all(param.is_valid() for param in self.parameters.values())
    
    def to_dict(self) -> Dict[str, Any]:
        """パラメータを辞書形式に変換"""
        return {name: param.value for name, param in self.parameters.items()}
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DockingParameters':
        """辞書からDockingParametersを作成"""
        instance = cls()
        for name, value in data.items():
            param_type = ParameterType.from_value(value)
            instance.add(DockingParameter(name=name, value=value, parameter_type=param_type))
        return instance
    
    def merge(self, other: 'DockingParameters') -> 'DockingParameters':
        """他のDockingParametersとマージした新しいインスタンスを作成
        
        他のパラメータが同じ名前のパラメータを持つ場合、それが優先されます。
        """
        result = DockingParameters()
        
        # 現在のパラメータをコピー
        for name, param in self.parameters.items():
            result.add(param)
        
        # 他のパラメータを追加（上書き）
        for name, param in other.parameters.items():
            result.add(param)
        
        return result