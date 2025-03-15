#!/usr/bin/env python3
"""
Ligandクラスのテスト

このモジュールはドッキング計算用のリガンド（準備された化合物）を表すLigandクラスをテストします。
"""

import pytest
from unittest.mock import MagicMock, patch
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.value_object.molecule_format import FormatType


class TestLigand:
    """Ligandクラスのテストケース"""
    
    @pytest.fixture
    def mock_compound(self):
        """テスト用のモック化合物"""
        compound = MagicMock(spec=Compound)
        compound.id = "C001"
        compound.name = "test_compound"
        compound.path = "test/compound.sdf"
        compound.metadata = {}
        compound.format = MagicMock()
        compound.format.type = FormatType.SDF
        # is_preparedプロパティをモック
        compound.is_prepared = False
        return compound
    
    def test_initialization(self, mock_compound):
        """初期化と基本プロパティのテスト"""
        # 最小限の引数での初期化
        ligand = Ligand(compound=mock_compound)
        
        assert ligand.compound is mock_compound
        assert ligand.metadata == {}
        assert ligand.name == mock_compound.id  # compoundのIDがデフォルトのname
        
        # 名前とメタデータを指定した初期化
        metadata = {"molecular_weight": 180.16}
        ligand = Ligand(
            compound=mock_compound,
            metadata=metadata,
            name="aspirin"
        )
        
        assert ligand.compound is mock_compound
        assert ligand.metadata is metadata
        assert ligand.name == "aspirin"
    
    def test_name_initialization_from_metadata(self, mock_compound):
        """メタデータからの名前初期化テスト"""
        # 化合物メタデータに名前がある場合
        mock_compound.metadata = {"name": "compound_in_metadata"}
        
        ligand = Ligand(compound=mock_compound)
        
        # メタデータの名前が優先される
        assert ligand.name == "compound_in_metadata"
        
        # 明示的に名前を指定した場合は、それが優先される
        ligand = Ligand(compound=mock_compound, name="explicit_name")
        
        assert ligand.name == "explicit_name"
    
    def test_is_prepared(self, mock_compound):
        """準備状態確認のテスト"""
        ligand = Ligand(compound=mock_compound)
        
        # 準備されていない状態（PDBQTではない）
        mock_compound.is_prepared = True
        mock_compound.format.type = FormatType.SDF
        assert ligand.is_prepared() is False
        
        # 準備されていない状態（is_preparedがFalse）
        mock_compound.is_prepared = False
        mock_compound.format.type = FormatType.PDBQT
        assert ligand.is_prepared() is False
        
        # 準備された状態
        mock_compound.is_prepared = True
        mock_compound.format.type = FormatType.PDBQT
        assert ligand.is_prepared() is True
    
    def test_get_path(self, mock_compound):
        """パス取得のテスト"""
        ligand = Ligand(compound=mock_compound)
        
        # compoundのパスが返される
        mock_compound.path = "test/path.sdf"
        assert ligand.get_path() == "test/path.sdf"
        
        # compoundのパスがNoneの場合
        mock_compound.path = None
        assert ligand.get_path() is None
    
    def test_get_id(self, mock_compound):
        """ID取得のテスト"""
        ligand = Ligand(compound=mock_compound)
        
        # compoundのIDが返される
        mock_compound.id = "test_id_123"
        assert ligand.get_id() == "test_id_123"
    
    def test_get_center_of_mass(self, mock_compound):
        """重心取得のテスト"""
        ligand = Ligand(compound=mock_compound)
        
        # compoundのget_center_of_massメソッドが呼ばれる
        mock_center = (1.0, 2.0, 3.0)
        mock_compound.get_center_of_mass.return_value = mock_center
        
        center = ligand.get_center_of_mass()
        
        assert center == mock_center
        mock_compound.get_center_of_mass.assert_called_once()
    
    def test_str_representation(self, mock_compound):
        """文字列表現のテスト"""
        ligand = Ligand(compound=mock_compound, name="test_ligand")
        
        # 準備状態がFalseの場合
        mock_compound.is_prepared = False
        mock_compound.format.type = FormatType.SDF
        
        str_rep = str(ligand)
        assert "test_ligand" in str_rep
        assert "prepared=False" in str_rep
        
        # 準備状態がTrueの場合
        mock_compound.is_prepared = True
        mock_compound.format.type = FormatType.PDBQT
        
        str_rep = str(ligand)
        assert "test_ligand" in str_rep
        assert "prepared=True" in str_rep