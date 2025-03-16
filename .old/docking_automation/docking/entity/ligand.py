from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Tuple

from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.value_object.molecule_format import FormatType


@dataclass
class Ligand:
    """ドッキング計算用のリガンド（準備された化合物）を表すエンティティ"""
    
    compound: Compound
    metadata: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None
    _is_prepared: bool = False
    
    def __post_init__(self) -> None:
        """初期化後の検証・処理"""
        # 名前が指定されていない場合、化合物のIDまたはメタデータから設定
        if not self.name:
            self.name = self.compound.metadata.get('name', self.compound.id)
    
    def is_prepared(self) -> bool:
        """リガンドが計算のために適切に準備されているかどうかを確認

        Returns:
            bool: リガンドが準備されているかどうか
        """
        # 化合物のフォーマットを確認
        if self.compound.format is None:
            return False
        
        # PDBQTフォーマットであることを確認
        format_type = getattr(self.compound.format, 'type', None)
        if format_type != FormatType.PDBQT:
            return False
            
        return self._is_prepared
    
    def set_prepared(self, is_prepared: bool) -> None:
        """リガンドの準備状態を設定する
        
        Args:
            is_prepared: 準備状態
        """
        self._is_prepared = is_prepared
    
    def get_path(self) -> Optional[str]:
        """リガンドファイルのパスを取得
        
        Returns:
            Optional[str]: ファイルパス
        """
        return self.compound.path
    
    def get_id(self) -> str:
        """リガンドの識別子を取得
        
        Returns:
            str: リガンドID
        """
        return self.compound.id
    
    def get_center_of_mass(self) -> Tuple[float, float, float]:
        """リガンドの重心を取得
        
        注意: 現在の実装では、デフォルト値を返します。
        実際の重心計算は別途実装が必要です。
        
        Returns:
            Tuple[float, float, float]: X, Y, Z座標
        """
        # TODO: 実際の重心計算を実装
        return (0.0, 0.0, 0.0)
    
    def __str__(self) -> str:
        """文字列表現
        
        Returns:
            str: インスタンスの文字列表現
        """
        return f"Ligand(name={self.name}, prepared={self.is_prepared()})"