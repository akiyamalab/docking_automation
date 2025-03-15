"""パラメータ型

このモジュールは、ドッキングパラメータの型を表す列挙型を提供します。
"""

from enum import Enum, auto
from typing import Any, Union, TypeVar, cast


T = TypeVar('T', int, float, bool, str)


class ParameterType(Enum):
    """ドッキングパラメータの型を表す列挙型
    
    ドッキングパラメータの型（整数、浮動小数点数、文字列、ブーリアンなど）を示します。
    """
    
    INTEGER = auto()  # 整数型
    FLOAT = auto()    # 浮動小数点数型
    STRING = auto()   # 文字列型
    BOOLEAN = auto()  # ブーリアン型
    
    @classmethod
    def from_value(cls, value: object) -> 'ParameterType':
        """値の型に基づいてパラメータ型を推定する
        
        Args:
            value: 型を推定する値
            
        Returns:
            推定されたパラメータ型
        """
        if isinstance(value, bool):
            return cls.BOOLEAN
        elif isinstance(value, int):
            return cls.INTEGER
        elif isinstance(value, float):
            return cls.FLOAT
        else:
            return cls.STRING
    
    def validate_value(self, value: object) -> bool:
        """値がこの型に適合するかどうかを検証する
        
        Args:
            value: 検証する値
            
        Returns:
            型に適合する場合はTrue、適合しない場合はFalse
        """
        if self is ParameterType.INTEGER:
            return isinstance(value, int) and not isinstance(value, bool)
        elif self is ParameterType.FLOAT:
            return isinstance(value, float)
        elif self is ParameterType.BOOLEAN:
            return isinstance(value, bool)
        elif self is ParameterType.STRING:
            return isinstance(value, str)
        return False
    
    def convert_value(self, value: Any) -> Union[int, float, bool, str]:
        """値をこの型に変換する
        
        可能であれば値をこの型に変換します。
        変換できない場合はValueErrorが発生します。
        
        Args:
            value: 変換する値
            
        Returns:
            変換された値
            
        Raises:
            ValueError: 値を変換できない場合
        """
        if self is ParameterType.INTEGER:
            return int(cast(Union[str, int, float], value))
        elif self is ParameterType.FLOAT:
            return float(cast(Union[str, int, float], value))
        elif self is ParameterType.BOOLEAN:
            if isinstance(value, str):
                value_str = value.lower()
                if value_str in ('true', 'yes', '1', 'on'):
                    return True
                if value_str in ('false', 'no', '0', 'off'):
                    return False
                raise ValueError(f"文字列 '{value_str}' をブーリアン値に変換できません")
            return bool(value)
        elif self is ParameterType.STRING:
            return str(value)
        
        raise ValueError(f"値 '{value}' を型 {self.name} に変換できません")