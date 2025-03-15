"""ドッキングパラメータ集合

このモジュールは、ドッキング計算のパラメータ集合を表現する値オブジェクトを提供します。
"""

from typing import Dict, List, Optional, Any, Iterator, Union

from docking_automation.docking.value_object.docking_parameter import DockingParameter


class DockingParameters:
    """ドッキング計算のパラメータ集合を表す値オブジェクト
    
    複数のドッキングパラメータをまとめて管理するためのコレクションクラスです。
    イテレータプロトコルをサポートし、含まれるパラメータを繰り返し処理できます。
    """
    
    def __init__(self, parameters: Optional[List[DockingParameter]] = None):
        """コンストラクタ
        
        Args:
            parameters: 初期パラメータのリスト（省略可）
        """
        self._parameters: Dict[str, DockingParameter] = {}
        
        if parameters:
            for param in parameters:
                self.add(param)
    
    def add(self, parameter: DockingParameter) -> None:
        """パラメータを追加する
        
        同じ名前のパラメータが既に存在する場合は上書きされます。
        
        Args:
            parameter: 追加するパラメータ
        """
        self._parameters[parameter.name] = parameter
    
    def get(self, name: str, default: Any = None) -> Union[DockingParameter, Any]:
        """名前を指定してパラメータを取得する
        
        Args:
            name: パラメータ名
            default: パラメータが存在しない場合のデフォルト値
            
        Returns:
            パラメータが存在する場合はそのパラメータ、存在しない場合はデフォルト値
        """
        return self._parameters.get(name, default)
    
    def remove(self, name: str) -> bool:
        """名前を指定してパラメータを削除する
        
        Args:
            name: 削除するパラメータの名前
            
        Returns:
            削除に成功した場合はTrue、パラメータが存在しなかった場合はFalse
        """
        if name in self._parameters:
            del self._parameters[name]
            return True
        return False
    
    def to_dict(self) -> Dict[str, Any]:
        """パラメータを辞書形式に変換する
        
        Returns:
            {パラメータ名: パラメータ値} 形式の辞書
        """
        return {name: param.value for name, param in self._parameters.items()}
    
    def __iter__(self) -> Iterator[DockingParameter]:
        """イテレータを返す
        
        Returns:
            パラメータを順に返すイテレータ
        """
        return iter(self._parameters.values())
    
    def __len__(self) -> int:
        """パラメータの数を返す
        
        Returns:
            含まれるパラメータの数
        """
        return len(self._parameters)
    
    def __bool__(self) -> bool:
        """空かどうかを判定する
        
        Returns:
            パラメータがある場合はTrue、ない場合はFalse
        """
        return bool(self._parameters)
    
    def __str__(self) -> str:
        """文字列表現を返す
        
        Returns:
            オブジェクトの文字列表現
        """
        return f"DockingParameters({', '.join(self._parameters.keys())})"