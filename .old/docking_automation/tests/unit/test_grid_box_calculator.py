"""グリッドボックス計算サービスのテスト"""

import pytest
import numpy as np
from typing import List, Tuple

from docking_automation.molecule.entity.protein import Protein
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.molecule.service.grid_box_calculator import GridBoxCalculator
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond


def create_test_structure() -> MoleculeStructure:
    """テスト用の分子構造を作成する
    
    簡単な立方体状の原子配置の構造を作成します
    
    Returns:
        MoleculeStructure: テスト用分子構造
    """
    atoms = [
        Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
        Atom(atom_id=2, element="C", x=10.0, y=0.0, z=0.0),
        Atom(atom_id=3, element="C", x=0.0, y=10.0, z=0.0),
        Atom(atom_id=4, element="C", x=10.0, y=10.0, z=0.0),
        Atom(atom_id=5, element="C", x=0.0, y=0.0, z=10.0),
        Atom(atom_id=6, element="C", x=10.0, y=0.0, z=10.0),
        Atom(atom_id=7, element="C", x=0.0, y=10.0, z=10.0),
        Atom(atom_id=8, element="C", x=10.0, y=10.0, z=10.0),
    ]
    
    bonds = [
        Bond(atom1_id=1, atom2_id=2, bond_type="single"),
        Bond(atom1_id=1, atom2_id=3, bond_type="single"),
        Bond(atom1_id=1, atom2_id=5, bond_type="single"),
        Bond(atom1_id=2, atom2_id=4, bond_type="single"),
        Bond(atom1_id=2, atom2_id=6, bond_type="single"),
        Bond(atom1_id=3, atom2_id=4, bond_type="single"),
        Bond(atom1_id=3, atom2_id=7, bond_type="single"),
        Bond(atom1_id=4, atom2_id=8, bond_type="single"),
        Bond(atom1_id=5, atom2_id=6, bond_type="single"),
        Bond(atom1_id=5, atom2_id=7, bond_type="single"),
        Bond(atom1_id=6, atom2_id=8, bond_type="single"),
        Bond(atom1_id=7, atom2_id=8, bond_type="single"),
    ]
    
    return MoleculeStructure(atoms=atoms, bonds=bonds)


def create_test_protein() -> Protein:
    """テスト用のタンパク質オブジェクトを作成する
    
    Returns:
        Protein: テスト用タンパク質
    """
    protein = Protein(
        id="test_protein",
        path="/path/to/test.pdb",
        chains=set(["A"])
    )
    
    # 構造を設定
    protein.structure = create_test_structure()
    
    # アクティブサイト情報を追加
    protein.metadata["active_site"] = {
        "residues": ["A:42", "A:45", "A:89"]
    }
    
    return protein


class TestGridBoxCalculator:
    """GridBoxCalculatorクラスのテスト"""
    
    def test_calculate_from_structure(self):
        """分子構造からグリッドボックスを計算するテスト"""
        # テスト用の構造を作成
        structure = create_test_structure()
        
        # デフォルトパラメータでグリッドボックスを計算
        grid_box = GridBoxCalculator.calculate_from_structure(structure)
        
        # 期待値：立方体の中心は (5, 5, 5)、サイズはパディング（デフォルト4.0Å）を含めて各方向 10 + 2*4 = 18Å
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(18.0)
        assert grid_box.size_y == pytest.approx(18.0)
        assert grid_box.size_z == pytest.approx(18.0)
        
        # パディングを変更してテスト
        padding = 2.0
        grid_box = GridBoxCalculator.calculate_from_structure(structure, padding=padding)
        
        # 期待値：立方体の中心は (5, 5, 5)、サイズはパディング2.0Åを含めて各方向 10 + 2*2 = 14Å
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(14.0)
        assert grid_box.size_y == pytest.approx(14.0)
        assert grid_box.size_z == pytest.approx(14.0)
        
        # 固定サイズで計算
        default_size = 20.0
        grid_box = GridBoxCalculator.calculate_from_structure(
            structure, 
            default_box_size=default_size
        )
        
        # 期待値：立方体の中心は (5, 5, 5)、サイズは各方向20Å
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(default_size)
        assert grid_box.size_y == pytest.approx(default_size)
        assert grid_box.size_z == pytest.approx(default_size)
    
    def test_calculate_from_protein(self):
        """タンパク質オブジェクトからグリッドボックスを計算するテスト"""
        # テスト用のタンパク質を作成
        protein = create_test_protein()
        
        # デフォルトパラメータでグリッドボックスを計算
        grid_box = GridBoxCalculator.calculate_from_protein(protein)
        
        # 期待値：デフォルトのボックスサイズ（20.0Å）が使用される
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(20.0)
        assert grid_box.size_y == pytest.approx(20.0)
        assert grid_box.size_z == pytest.approx(20.0)
        
        # 固定サイズなしでグリッドボックスを計算
        grid_box = GridBoxCalculator.calculate_from_protein(
            protein, 
            default_box_size=None
        )
        
        # 期待値：実際のサイズ + パディング（4.0Å）
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(18.0)
        assert grid_box.size_y == pytest.approx(18.0)
        assert grid_box.size_z == pytest.approx(18.0)
    
    def test_protein_calculate_grid_box(self):
        """Proteinクラスのcalculate_grid_boxメソッドのテスト"""
        # テスト用のタンパク質を作成
        protein = create_test_protein()
        
        # デフォルトパラメータでグリッドボックスを計算
        grid_box = protein.calculate_grid_box()
        
        # 期待値：デフォルトのボックスサイズ（20.0Å）が使用される
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        assert grid_box.size_x == pytest.approx(20.0)
        assert grid_box.size_y == pytest.approx(20.0)
        assert grid_box.size_z == pytest.approx(20.0)
        
        # アクティブサイトのみを使用してグリッドボックスを計算
        # 注：現在の実装では仮実装のため、結果は同じになる
        grid_box = protein.calculate_grid_box(active_site_only=True)
        
        # アクティブサイトから計算する処理はモック実装なので、
        # 現状では全体と同じ結果になるが、将来的には異なる結果になるはず
        assert grid_box.center_x == pytest.approx(5.0)
        assert grid_box.center_y == pytest.approx(5.0)
        assert grid_box.center_z == pytest.approx(5.0)
        
    def test_get_bounding_box(self):
        """MoleculeStructureのget_bounding_boxメソッドのテスト"""
        # テスト用の構造を作成
        structure = create_test_structure()
        
        # 境界ボックスを取得
        min_x, min_y, min_z, max_x, max_y, max_z = structure.get_bounding_box()
        
        # 期待値：原子座標の最小/最大値
        assert min_x == pytest.approx(0.0)
        assert min_y == pytest.approx(0.0)
        assert min_z == pytest.approx(0.0)
        assert max_x == pytest.approx(10.0)
        assert max_y == pytest.approx(10.0)
        assert max_z == pytest.approx(10.0)
        
    def test_get_center_of_mass(self):
        """MoleculeStructureのget_center_of_massメソッドのテスト"""
        # テスト用の構造を作成
        structure = create_test_structure()
        
        # 重心を取得
        center_x, center_y, center_z = structure.get_center_of_mass()
        
        # 期待値：立方体の中心（5.0, 5.0, 5.0）
        assert center_x == pytest.approx(5.0)
        assert center_y == pytest.approx(5.0)
        assert center_z == pytest.approx(5.0)