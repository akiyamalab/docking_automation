from dataclasses import dataclass, field
from typing import Dict, Any, Optional, Set

from docking_automation.domain.molecule.entity.protein import Protein
from docking_automation.domain.molecule.value_object.molecule_format import FormatType
from docking_automation.domain.docking.value_object.grid_box import GridBox


@dataclass
class Receptor:
    """ドッキング計算用のレセプター（準備されたタンパク質）を表すエンティティ"""
    
    protein: Protein
    active_site: Optional[GridBox] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    name: Optional[str] = None
    
    def __post_init__(self) -> None:
        """初期化後の検証・処理"""
        # 名前が指定されていない場合、タンパク質のIDまたはメタデータから設定
        if not self.name:
            self.name = self.protein.metadata.get('name', self.protein.id)
        
        # アクティブサイトが指定されていない場合で、タンパク質にアクティブサイト残基が
        # 設定されている場合は、それらの残基からグリッドボックスを計算
        if not self.active_site and self.protein.active_site_residues:
            center = self.protein.get_active_site_center()
            
            # デフォルトサイズのグリッドボックスを作成
            self.active_site = GridBox(
                center_x=center[0],
                center_y=center[1],
                center_z=center[2],
                size_x=20.0,
                size_y=20.0,
                size_z=20.0
            )
    
    def is_prepared(self) -> bool:
        """レセプターが計算のために適切に準備されているかどうかを確認"""
        # タンパク質が準備されていて、PDBQTフォーマットであることを確認
        return (
            self.protein.is_prepared and 
            self.protein.format.type == FormatType.PDBQT
        )
    
    def get_path(self) -> Optional[str]:
        """レセプターファイルのパスを取得"""
        return self.protein.path
    
    def get_id(self) -> str:
        """レセプターの識別子を取得"""
        return self.protein.id
    
    def get_chains(self) -> Set[str]:
        """レセプターのチェーンIDを取得"""
        return self.protein.chains
    
    def get_suggested_grid_box(self) -> GridBox:
        """提案されるグリッドボックスを取得
        
        アクティブサイトが設定されている場合はそれを返し、
        そうでない場合はタンパク質の中心に標準サイズのグリッドボックスを返す
        """
        if self.active_site:
            return self.active_site
        
        # タンパク質の中心に標準サイズのグリッドボックスを配置
        center = self.protein.get_center_of_mass()
        return GridBox(
            center_x=center[0],
            center_y=center[1],
            center_z=center[2],
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
    
    def has_active_site_defined(self) -> bool:
        """アクティブサイトが明示的に定義されているかどうかを確認"""
        return self.active_site is not None or bool(self.protein.active_site_residues)
    
    def with_active_site(self, grid_box: GridBox) -> 'Receptor':
        """指定されたアクティブサイトで新しいレセプターを作成"""
        return Receptor(
            protein=self.protein,
            active_site=grid_box,
            metadata=self.metadata.copy(),
            name=self.name
        )
    
    def __str__(self) -> str:
        """文字列表現"""
        return f"Receptor(name={self.name}, prepared={self.is_prepared()}, active_site_defined={self.has_active_site_defined()})"