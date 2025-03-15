from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any

from docking_automation.domain.molecule.entity.molecule import Molecule
from docking_automation.domain.molecule.value_object.molecule_structure import MoleculeStructure
from docking_automation.domain.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.domain.molecule.value_object.molecule_property import MoleculeProperty


@dataclass
class Compound(Molecule):
    """化合物（リガンド）を表すエンティティ"""
    
    # Molecule基底クラスから継承するフィールド:
    # id: str
    # structure: MoleculeStructure
    # format: MoleculeFormat
    # properties: List[MoleculeProperty]
    # path: Optional[str]
    # metadata: Dict[str, Any]
    
    is_prepared: bool = False  # ドッキング計算のための準備が完了しているかどうか
    preparation_method: Optional[str] = None  # 準備に使用したメソッド（例: "gasteiger", "am1bcc"）
    
    def prepare(self) -> 'Compound':
        """ドッキング計算のための準備を行う
        
        このメソッドは、インフラストラクチャ層の実際の実装で具体的な準備処理を行います。
        ドメイン層では抽象的なインターフェースのみを定義します。
        """
        if self.is_prepared:
            return self
        
        # 具体的な準備処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして新しいインスタンスを返す
        prepared_compound = Compound(
            id=self.id,
            structure=self.structure,
            format=self.format,
            properties=self.properties.copy(),
            path=self.path,
            metadata=self.metadata.copy(),
            is_prepared=True,
            preparation_method="default"
        )
        
        return prepared_compound
    
    def convert_to(self, target_format: FormatType) -> 'Compound':
        """指定された形式に変換する
        
        Args:
            target_format: 変換先の形式
            
        Returns:
            変換された新しいCompoundインスタンス
        """
        # 同じ形式の場合は自分自身を返す
        if self.format.type == target_format:
            return self
        
        # 変換後の形式を表すMoleculeFormatオブジェクト
        new_format = MoleculeFormat(type=target_format)
        
        # 具体的な変換処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして新しいインスタンスを返す
        converted_compound = Compound(
            id=self.id,
            structure=self.structure,
            format=new_format,
            properties=self.properties.copy(),
            path=None,  # 変換後のファイルパスはまだ存在しない
            metadata=self.metadata.copy(),
            is_prepared=self.is_prepared,
            preparation_method=self.preparation_method
        )
        
        return converted_compound
    
    def to_ligand(self) -> 'Compound':
        """ドッキング計算用のリガンドとして準備された化合物を返す
        
        PDBQTフォーマットへの変換と必要な準備処理を行います。
        
        Returns:
            ドッキング計算用に準備された化合物インスタンス
        """
        # まだ準備されていない場合は準備する
        if not self.is_prepared:
            prepared = self.prepare()
        else:
            prepared = self
        
        # PDBQTに変換
        if prepared.format.type != FormatType.PDBQT:
            prepared = prepared.convert_to(FormatType.PDBQT)
        
        return prepared
    
    def __str__(self) -> str:
        """文字列表現"""
        name = self.metadata.get('name', self.id)
        return f"Compound(id={self.id}, name={name}, format={self.format.type.name}, prepared={self.is_prepared})"