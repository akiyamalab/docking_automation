from dataclasses import dataclass, field
from typing import List, Optional, Dict, Any, Set, Tuple

from docking_automation.domain.molecule.entity.molecule import Molecule
from docking_automation.domain.molecule.value_object.molecule_structure import MoleculeStructure
from docking_automation.domain.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.domain.molecule.value_object.molecule_property import MoleculeProperty


@dataclass
class Protein(Molecule):
    """タンパク質（レセプター）を表すエンティティ"""
    
    # Molecule基底クラスから継承するフィールド:
    # id: str
    # structure: MoleculeStructure
    # format: MoleculeFormat
    # properties: List[MoleculeProperty]
    # path: Optional[str]
    # metadata: Dict[str, Any]
    
    is_prepared: bool = False  # ドッキング計算のための準備が完了しているかどうか
    preparation_method: Optional[str] = None  # 準備に使用したメソッド
    chains: Set[str] = field(default_factory=set)  # タンパク質のチェーンID
    has_water: bool = True  # 水分子を含むかどうか
    has_hydrogens: bool = False  # 水素原子を含むかどうか
    active_site_residues: Set[str] = field(default_factory=set)  # 活性部位残基のID
    
    def prepare(self) -> 'Protein':
        """ドッキング計算のための準備を行う
        
        以下の処理を行います：
        1. 水分子の除去（オプション）
        2. 水素原子の追加
        3. 電荷の計算と割り当て
        
        Returns:
            準備されたタンパク質インスタンス
        """
        if self.is_prepared:
            return self
        
        # 具体的な準備処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして新しいインスタンスを返す
        prepared_protein = Protein(
            id=self.id,
            structure=self.structure,
            format=self.format,
            properties=self.properties.copy(),
            path=self.path,
            metadata=self.metadata.copy(),
            is_prepared=True,
            preparation_method="default",
            chains=self.chains.copy(),
            has_water=False,  # 水分子は除去されたと仮定
            has_hydrogens=True,  # 水素原子は追加されたと仮定
            active_site_residues=self.active_site_residues.copy()
        )
        
        return prepared_protein
    
    def convert_to(self, target_format: FormatType) -> 'Protein':
        """指定された形式に変換する
        
        Args:
            target_format: 変換先の形式
            
        Returns:
            変換された新しいProteinインスタンス
        """
        # 同じ形式の場合は自分自身を返す
        if self.format.type == target_format:
            return self
        
        # 変換後の形式を表すMoleculeFormatオブジェクト
        new_format = MoleculeFormat(type=target_format)
        
        # 具体的な変換処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして新しいインスタンスを返す
        converted_protein = Protein(
            id=self.id,
            structure=self.structure,
            format=new_format,
            properties=self.properties.copy(),
            path=None,  # 変換後のファイルパスはまだ存在しない
            metadata=self.metadata.copy(),
            is_prepared=self.is_prepared,
            preparation_method=self.preparation_method,
            chains=self.chains.copy(),
            has_water=self.has_water,
            has_hydrogens=self.has_hydrogens,
            active_site_residues=self.active_site_residues.copy()
        )
        
        return converted_protein
    
    def remove_water(self) -> 'Protein':
        """水分子を除去した新しいタンパク質インスタンスを返す"""
        if not self.has_water:
            return self
        
        # 水分子除去の具体的な処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして水分子フラグを変更した新しいインスタンスを返す
        protein_without_water = Protein(
            id=self.id,
            structure=self.structure,  # 実際には水分子を除去した構造になるはず
            format=self.format,
            properties=self.properties.copy(),
            path=self.path,
            metadata=self.metadata.copy(),
            is_prepared=self.is_prepared,
            preparation_method=self.preparation_method,
            chains=self.chains.copy(),
            has_water=False,
            has_hydrogens=self.has_hydrogens,
            active_site_residues=self.active_site_residues.copy()
        )
        
        return protein_without_water
    
    def add_hydrogens(self) -> 'Protein':
        """水素原子を追加した新しいタンパク質インスタンスを返す"""
        if self.has_hydrogens:
            return self
        
        # 水素追加の具体的な処理はインフラストラクチャ層で行われるため、
        # ここではプレースホルダーとして水素フラグを変更した新しいインスタンスを返す
        protein_with_hydrogens = Protein(
            id=self.id,
            structure=self.structure,  # 実際には水素を追加した構造になるはず
            format=self.format,
            properties=self.properties.copy(),
            path=self.path,
            metadata=self.metadata.copy(),
            is_prepared=self.is_prepared,
            preparation_method=self.preparation_method,
            chains=self.chains.copy(),
            has_water=self.has_water,
            has_hydrogens=True,
            active_site_residues=self.active_site_residues.copy()
        )
        
        return protein_with_hydrogens
    
    def to_receptor(self) -> 'Protein':
        """ドッキング計算用のレセプターとして準備されたタンパク質を返す
        
        PDBQTフォーマットへの変換と必要な準備処理を行います。
        
        Returns:
            ドッキング計算用に準備されたタンパク質インスタンス
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
    
    def get_active_site_center(self) -> Tuple[float, float, float]:
        """活性部位の中心座標を計算する
        
        活性部位が定義されていない場合はタンパク質全体の重心を返す
        
        Returns:
            活性部位の中心座標 (x, y, z)
        """
        # 活性部位が定義されていない場合はタンパク質全体の重心を返す
        if not self.active_site_residues:
            return self.get_center_of_mass()
        
        # 実際の実装では、活性部位残基に属する原子の座標から中心を計算する
        # ここではプレースホルダーとして重心を返す
        return self.get_center_of_mass()
    
    def __str__(self) -> str:
        """文字列表現"""
        name = self.metadata.get('name', self.id)
        chain_str = ','.join(self.chains) if self.chains else 'None'
        return f"Protein(id={self.id}, name={name}, chains={chain_str}, format={self.format.type.name}, prepared={self.is_prepared})"