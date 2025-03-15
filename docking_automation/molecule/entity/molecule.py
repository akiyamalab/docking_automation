from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Tuple
import uuid

from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure
from docking_automation.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.molecule.value_object.molecule_property import MoleculeProperty, MoleculeProperties


@dataclass
class Molecule(ABC):
    """分子の抽象基底クラス"""
    id: str  # 分子の一意識別子
    structure: MoleculeStructure
    format: MoleculeFormat
    properties: List[MoleculeProperty] = field(default_factory=list)
    path: Optional[str] = None  # 分子のファイルパス（存在する場合）
    metadata: Dict[str, Any] = field(default_factory=dict)  # 追加のメタデータ
    
    def __post_init__(self) -> None:
        """初期化後の処理"""
        # IDが指定されていない場合はUUIDを生成
        if not self.id:
            self.id = uuid.uuid4().hex
    
    @abstractmethod
    def prepare(self) -> 'Molecule':
        """ドッキング計算のための準備を行う（抽象メソッド）"""
        pass
    
    def convert_to(self, target_format: FormatType) -> 'Molecule':
        """指定された形式に変換した新しい分子インスタンスを返す"""
        # 同じ形式の場合は自分自身を返す
        if self.format.type == target_format:
            return self
        
        # このメソッドはサブクラスでオーバーライドされる可能性があるため、
        # デフォルトではNotImplementedErrorをスロー
        raise NotImplementedError(
            f"Conversion from {self.format.type} to {target_format} is not implemented for {type(self).__name__}"
        )
    
    def validate(self) -> bool:
        """分子の構造と形式が有効かどうかを検証"""
        # 構造の検証
        if not self.structure.validate():
            return False
        
        # 形式の検証
        if not self.format.validate():
            return False
        
        # プロパティの検証
        for prop in self.properties:
            if not prop.is_valid():
                return False
        
        return True
    
    def get_property(self, name: str) -> Optional[MoleculeProperty]:
        """指定した名前のプロパティを取得"""
        for prop in self.properties:
            if prop.name == name:
                return prop
        return None
    
    def get_all_properties(self) -> MoleculeProperties:
        """全てのプロパティをMoleculePropertiesオブジェクトとして取得"""
        return MoleculeProperties(properties=self.properties)
    
    def add_property(self, property: MoleculeProperty) -> None:
        """プロパティを追加（既存の同名プロパティは上書き）"""
        # 既存の同名プロパティを削除
        self.properties = [p for p in self.properties if p.name != property.name]
        # 新しいプロパティを追加
        self.properties.append(property)
    
    def get_center_of_mass(self) -> Tuple[float, float, float]:
        """分子の重心を取得"""
        return self.structure.get_center_of_mass()
    
    def get_file_name(self) -> str:
        """分子のファイル名を取得（存在しない場合はIDを使用）"""
        if self.path:
            return self.path.split('/')[-1]
        return f"{self.id}.{self.format.type.get_extension()}"