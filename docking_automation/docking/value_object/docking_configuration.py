from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Union, cast

from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.docking_parameter import DockingParameters


@dataclass(frozen=True)
class DockingConfiguration:
    """ドッキング計算の設定を表す値オブジェクト
    
    ドッキング計算に必要な設定情報（グリッドボックス、パラメータなど）を表します。
    """
    grid_box: GridBox
    parameters: DockingParameters
    name: Optional[str] = None
    description: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def validate(self) -> bool:
        """設定が有効かどうかを検証"""
        # グリッドボックスの検証
        try:
            # GridBoxコンストラクタ内で検証が行われる
            GridBox(**self.grid_box.to_dict())
        except ValueError:
            return False
        
        # パラメータの検証
        if not self.parameters.validate():
            return False
        
        return True
    
    def with_parameters(self, parameters: DockingParameters) -> 'DockingConfiguration':
        """指定されたパラメータを適用した新しい設定オブジェクトを作成"""
        merged_params = self.parameters.merge(parameters)
        
        return DockingConfiguration(
            grid_box=self.grid_box,
            parameters=merged_params,
            name=self.name,
            description=self.description,
            metadata=self.metadata.copy()
        )
    
    def with_grid_box(self, grid_box: GridBox) -> 'DockingConfiguration':
        """指定されたグリッドボックスを適用した新しい設定オブジェクトを作成"""
        return DockingConfiguration(
            grid_box=grid_box,
            parameters=self.parameters,
            name=self.name,
            description=self.description,
            metadata=self.metadata.copy()
        )
    
    def get_parameter(self, name: str, default: Any = None) -> Any:
        """パラメータの値を取得"""
        return self.parameters.get_value(name, default)
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        result = {
            'grid_box': self.grid_box.to_dict(),
            'parameters': self.parameters.to_dict(),
        }
        
        if self.name:
            result['name'] = cast(Any, self.name)
            
        if self.description:
            result['description'] = cast(Any, self.description)
            
        if self.metadata:
            result['metadata'] = self.metadata
            
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'DockingConfiguration':
        """辞書からDockingConfigurationを作成"""
        grid_box = GridBox.from_dict(data['grid_box'])
        parameters = DockingParameters.from_dict(data['parameters'])
        
        return cls(
            grid_box=grid_box,
            parameters=parameters,
            name=data.get('name'),
            description=data.get('description'),
            metadata=data.get('metadata', {})
        )
    
    @classmethod
    def create_default(cls, center_x: float, center_y: float, center_z: float) -> 'DockingConfiguration':
        """デフォルト設定を作成
        
        Args:
            center_x: グリッドボックスのX中心座標
            center_y: グリッドボックスのY中心座標
            center_z: グリッドボックスのZ中心座標
            
        Returns:
            デフォルト設定のDockingConfiguration
        """
        grid_box = GridBox(
            center_x=center_x,
            center_y=center_y,
            center_z=center_z,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
        
        # AutoDock Vinaのデフォルトパラメータ
        parameters = DockingParameters.from_dict({
            'exhaustiveness': 8,
            'num_modes': 9,
            'energy_range': 3,
            'seed': 0,
            'cpu': 1
        })
        
        return cls(
            grid_box=grid_box,
            parameters=parameters,
            name="Default Configuration"
        )