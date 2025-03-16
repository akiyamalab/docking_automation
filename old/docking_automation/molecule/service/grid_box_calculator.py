"""グリッドボックス計算サービス

このモジュールは、分子構造からドッキング計算のためのグリッドボックスを計算するためのサービスクラスを提供します。
"""

from typing import Optional, Tuple, Union, List

from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.entity.molecule import Molecule
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure


class GridBoxCalculator:
    """分子構造からグリッドボックスを計算するためのサービスクラス
    
    タンパク質構造やリガンド構造から、ドッキング計算に適したグリッドボックスの
    中心座標とサイズを計算します。
    """
    
    @staticmethod
    def calculate_from_structure(
        structure: MoleculeStructure,
        padding: float = 4.0,
        default_box_size: Optional[float] = None
    ) -> GridBox:
        """分子構造からグリッドボックスを計算する
        
        Args:
            structure: 分子構造オブジェクト
            padding: 境界ボックスに追加するパディング（Å）
            default_box_size: ボックスサイズを固定値にする場合の値（Å）
        
        Returns:
            計算されたグリッドボックス
        
        Examples:
            >>> from docking_automation.molecule.entity.protein import Protein
            >>> protein = Protein.from_file("path/to/protein.pdb")
            >>> grid_box = GridBoxCalculator.calculate_from_structure(protein.structure)
        """
        # 分子の中心を計算
        center = structure.get_center_of_mass()
        
        if default_box_size is not None:
            # 固定サイズのボックスを作成
            size = (default_box_size, default_box_size, default_box_size)
        else:
            # 分子の境界ボックスを取得
            min_x, min_y, min_z, max_x, max_y, max_z = structure.get_bounding_box()
            
            # 境界ボックスのサイズを計算（パディング付き）
            size_x = max_x - min_x + 2 * padding
            size_y = max_y - min_y + 2 * padding
            size_z = max_z - min_z + 2 * padding
            size = (size_x, size_y, size_z)
        
        return GridBox(
            center_x=center[0],
            center_y=center[1],
            center_z=center[2],
            size_x=size[0],
            size_y=size[1],
            size_z=size[2]
        )
    
    @staticmethod
    def calculate_from_molecule(
        molecule: Molecule,
        padding: float = 4.0,
        default_box_size: Optional[float] = None
    ) -> GridBox:
        """分子オブジェクトからグリッドボックスを計算する
        
        Args:
            molecule: 分子オブジェクト
            padding: 境界ボックスに追加するパディング（Å）
            default_box_size: ボックスサイズを固定値にする場合の値（Å）
        
        Returns:
            計算されたグリッドボックス
        """
        return GridBoxCalculator.calculate_from_structure(
            molecule.structure,
            padding=padding,
            default_box_size=default_box_size
        )
    
    @staticmethod
    def calculate_from_protein(
        protein: Protein,
        active_site_residues: Optional[List[str]] = None,
        padding: float = 4.0,
        default_box_size: Optional[float] = 20.0
    ) -> GridBox:
        """タンパク質オブジェクトからグリッドボックスを計算する
        
        Args:
            protein: タンパク質オブジェクト
            active_site_residues: 活性部位の残基ID（指定がある場合はそれらを中心にボックスを作成）
            padding: 境界ボックスに追加するパディング（Å）
            default_box_size: デフォルトのボックスサイズ（Å）
        
        Returns:
            計算されたグリッドボックス
            
        Notes:
            active_site_residuesが指定されている場合、それらの残基を中心にボックスを作成します。
            指定がない場合は、タンパク質全体から計算されます。その場合、default_box_sizeが
            指定されていればその値が使用されます。
        """
        # タンパク質の構造がない場合はエラー
        if not hasattr(protein, 'structure') or protein.structure is None:
            raise ValueError("タンパク質の構造情報がありません")
        
        # アクティブサイト残基が指定されている場合、それらの原子のみを抽出
        # （このままでは動作しないため、Proteinクラスの改修が必要）
        if active_site_residues and hasattr(protein, 'get_residue_atoms'):
            # この部分は仮実装で、Proteinクラスの拡張が必要
            active_site_atoms = protein.get_residue_atoms(active_site_residues)
            if active_site_atoms:
                # アクティブサイト原子から構造を作成して計算
                # この部分は実際のProteinクラスの実装に合わせて変更が必要
                return GridBoxCalculator.calculate_from_structure(
                    active_site_atoms,
                    padding=padding,
                    default_box_size=None  # アクティブサイトの場合は実際のサイズを使用
                )
        
        # 通常のケース：タンパク質全体から計算
        return GridBoxCalculator.calculate_from_structure(
            protein.structure,
            padding=padding,
            default_box_size=default_box_size
        )
    
    @staticmethod
    def calculate_from_ligand_protein(
        ligand: Molecule,
        protein: Protein,
        padding: float = 4.0
    ) -> GridBox:
        """リガンドとタンパク質からグリッドボックスを計算する
        
        リガンドの位置を中心に、適切なサイズのグリッドボックスを計算します。
        結晶構造複合体からドッキング領域を定義する場合に使用します。
        
        Args:
            ligand: リガンド分子
            protein: タンパク質
            padding: リガンド周囲に追加するパディング（Å）
        
        Returns:
            リガンドを中心としたグリッドボックス
        """
        # リガンドからグリッドボックスを計算
        ligand_box = GridBoxCalculator.calculate_from_molecule(
            ligand, 
            padding=padding,
            default_box_size=None  # リガンドの実際のサイズを使用
        )
        
        # ボックスが小さすぎる場合は最小サイズを保証
        min_size = 10.0  # 最小ボックスサイズ（Å）
        size_x, size_y, size_z = ligand_box.get_size()
        
        if size_x < min_size:
            size_x = min_size
        if size_y < min_size:
            size_y = min_size
        if size_z < min_size:
            size_z = min_size
        
        # 調整されたサイズで新しいGridBoxを作成
        return GridBox(
            center_x=ligand_box.center_x,
            center_y=ligand_box.center_y,
            center_z=ligand_box.center_z,
            size_x=size_x,
            size_y=size_y,
            size_z=size_z
        )