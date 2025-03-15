from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Tuple

from docking_automation.domain.molecule.entity.compound import Compound
from docking_automation.domain.molecule.value_object.molecule_format import FormatType


@dataclass
class Ligand:
    """ドッキング計算用のリガンド（準備された化合物）を表すエンティティ"""
    
    compound: Compound
    metadata: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None
    
    def __post_init__(self) -> None:
        """初期化後の検証・処理"""
        # 名前が指定されていない場合、化合物のIDまたはメタデータから設定
        if not self.name:
            self.name = self.compound.metadata.get('name', self.compound.id)
    
    def is_prepared(self) -> bool:
        """リガンドが計算のために適切に準備されているかどうかを確認"""
        # 化合物が準備されていて、PDBQTフォーマットであることを確認
        return (
            self.compound.is_prepared and 
            self.compound.format.type == FormatType.PDBQT
        )
    
    def get_path(self) -> Optional[str]:
        """リガンドファイルのパスを取得"""
        return self.compound.path
    
    def get_id(self) -> str:
        """リガンドの識別子を取得"""
        return self.compound.id
    
    def get_center_of_mass(self) -> Tuple[float, float, float]:
        """リガンドの重心を取得"""
        return self.compound.get_center_of_mass()
    
    def __str__(self) -> str:
        """文字列表現"""
        return f"Ligand(name={self.name}, prepared={self.is_prepared()})"