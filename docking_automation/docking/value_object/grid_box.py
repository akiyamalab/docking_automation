from dataclasses import dataclass
from typing import Optional, Tuple, Dict, Any


@dataclass(frozen=True)
class GridBox:
    """ドッキング計算の探索空間を表す値オブジェクト
    
    ドッキング計算において、タンパク質のどの領域に対して化合物の結合ポーズを探索するかを
    定義する3次元の箱（ボックス）を表します。
    """
    
    # 中心座標
    center_x: float
    center_y: float
    center_z: float
    
    # ボックスのサイズ
    size_x: float
    size_y: float
    size_z: float
    
    # 追加情報（オプション）
    spacing: float = 1.0  # グリッドポイント間の距離（Å）
    
    def __post_init__(self) -> None:
        """初期化後の検証"""
        if self.size_x <= 0 or self.size_y <= 0 or self.size_z <= 0:
            raise ValueError("Grid box sizes must be positive")
        
        if self.spacing <= 0:
            raise ValueError("Grid spacing must be positive")
    
    def get_center(self) -> Tuple[float, float, float]:
        """ボックスの中心座標を取得"""
        return (self.center_x, self.center_y, self.center_z)
    
    def get_size(self) -> Tuple[float, float, float]:
        """ボックスのサイズを取得"""
        return (self.size_x, self.size_y, self.size_z)
    
    def get_corner_min(self) -> Tuple[float, float, float]:
        """ボックスの最小座標を取得"""
        half_x = self.size_x / 2
        half_y = self.size_y / 2
        half_z = self.size_z / 2
        
        return (
            self.center_x - half_x,
            self.center_y - half_y,
            self.center_z - half_z
        )
    
    def get_corner_max(self) -> Tuple[float, float, float]:
        """ボックスの最大座標を取得"""
        half_x = self.size_x / 2
        half_y = self.size_y / 2
        half_z = self.size_z / 2
        
        return (
            self.center_x + half_x,
            self.center_y + half_y,
            self.center_z + half_z
        )
    
    def contains_point(self, x: float, y: float, z: float) -> bool:
        """指定された点がボックス内に含まれるかどうかを判定"""
        min_x, min_y, min_z = self.get_corner_min()
        max_x, max_y, max_z = self.get_corner_max()
        
        return (
            min_x <= x <= max_x and
            min_y <= y <= max_y and
            min_z <= z <= max_z
        )
    
    def get_volume(self) -> float:
        """ボックスの体積を計算"""
        return self.size_x * self.size_y * self.size_z
    
    def with_padding(self, padding: float) -> 'GridBox':
        """指定されたパディングを加えた新しいグリッドボックスを作成"""
        if padding < 0:
            raise ValueError("Padding must be non-negative")
        
        return GridBox(
            center_x=self.center_x,
            center_y=self.center_y,
            center_z=self.center_z,
            size_x=self.size_x + 2 * padding,
            size_y=self.size_y + 2 * padding,
            size_z=self.size_z + 2 * padding,
            spacing=self.spacing
        )
    
    def to_dict(self) -> Dict[str, float]:
        """辞書形式に変換"""
        return {
            'center_x': self.center_x,
            'center_y': self.center_y,
            'center_z': self.center_z,
            'size_x': self.size_x,
            'size_y': self.size_y,
            'size_z': self.size_z,
            'spacing': self.spacing
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'GridBox':
        """辞書からGridBoxを作成"""
        return cls(
            center_x=float(data['center_x']),
            center_y=float(data['center_y']),
            center_z=float(data['center_z']),
            size_x=float(data['size_x']),
            size_y=float(data['size_y']),
            size_z=float(data['size_z']),
            spacing=float(data.get('spacing', 1.0))
        )