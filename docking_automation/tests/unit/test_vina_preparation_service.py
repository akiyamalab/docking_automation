#!/usr/bin/env python3
"""
VinaPreparationServiceクラスのテスト

このモジュールはAutoDock Vina用の分子準備サービスをテストします。
VinaPreparationServiceクラスは、AutoDock Vinaでのドッキング計算に必要な分子の準備を行います。
"""

import pytest
from unittest.mock import MagicMock
from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService
from docking_automation.molecule.service.molecule_preparation_service import MoleculeProperty
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein


class TestVinaPreparationService:
    """VinaPreparationServiceクラスのテストケース"""
    
    @pytest.fixture
    def service(self):
        """テスト用のVinaPreparationServiceインスタンス"""
        return VinaPreparationService()
    
    @pytest.fixture
    def mock_compound(self):
        """テスト用のモック化合物"""
        compound = MagicMock(spec=Compound)
        compound.id = "C001"
        compound.name = "test_compound"
        compound.path = "test/compound.sdf"
        compound.metadata = {}
        return compound
    
    @pytest.fixture
    def mock_protein(self):
        """テスト用のモックタンパク質"""
        protein = MagicMock(spec=Protein)
        protein.id = "P001"
        protein.name = "test_protein"
        protein.path = "test/protein.pdb"
        protein.chains = {"A", "B"}
        protein.metadata = {}
        return protein
    
    def test_prepare_ligand(self, service, mock_compound):
        """prepare_ligand()メソッドのテスト"""
        # デフォルトメソッドでの準備
        prepared_compound = service.prepare_ligand(mock_compound)
        
        # 同じインスタンスが返されることを確認
        assert prepared_compound is mock_compound
        
        # メタデータが追加されていることを確認
        mock_compound.add_metadata.assert_called_once()
        args, _ = mock_compound.add_metadata.call_args
        metadata = args[0]
        assert "preparation_tool" in metadata
        assert metadata["preparation_tool"] == "autodock_vina"
        assert metadata["preparation_method"] == "default"
        
        # 特定のメソッドを指定した場合
        mock_compound.add_metadata.reset_mock()
        
        prepared_compound = service.prepare_ligand(mock_compound, method="minimize")
        
        mock_compound.add_metadata.assert_called_once()
        args, _ = mock_compound.add_metadata.call_args
        metadata = args[0]
        assert metadata["preparation_method"] == "minimize"
    
    def test_prepare_receptor(self, service, mock_protein):
        """prepare_receptor()メソッドのテスト"""
        # デフォルトメソッドでの準備（水分子除去あり）
        with pytest.monkeypatch.context() as m:
            # remove_waterメソッドのモック化
            m.setattr(service, "remove_water", lambda p: p)
            
            prepared_protein = service.prepare_receptor(mock_protein)
            
            # 同じインスタンスが返されることを確認
            assert prepared_protein is mock_protein
            
            # メタデータが追加されていることを確認
            mock_protein.add_metadata.assert_called_once()
            args, _ = mock_protein.add_metadata.call_args
            metadata = args[0]
            assert "preparation_tool" in metadata
            assert metadata["preparation_tool"] == "autodock_vina"
            assert metadata["preparation_method"] == "default"
        
        # 水分子を保持するメソッドでの準備
        mock_protein.add_metadata.reset_mock()
        
        # remove_waterメソッドがスパイ
        remove_water_spy = MagicMock(return_value=mock_protein)
        service.remove_water = remove_water_spy
        
        prepared_protein = service.prepare_receptor(mock_protein)
        
        # remove_waterメソッドが呼ばれたことを確認
        remove_water_spy.assert_called_once_with(mock_protein)
        
        # 水分子を保持するメソッドの場合
        remove_water_spy.reset_mock()
        mock_protein.add_metadata.reset_mock()
        
        prepared_protein = service.prepare_receptor(mock_protein, method="keepwater")
        
        # remove_waterメソッドが呼ばれていないことを確認
        remove_water_spy.assert_not_called()
        
        # メタデータが追加されていることを確認
        mock_protein.add_metadata.assert_called_once()
        args, _ = mock_protein.add_metadata.call_args
        metadata = args[0]
        assert metadata["preparation_method"] == "keepwater"
    
    def test_calculate_properties(self, service, mock_compound):
        """calculate_properties()メソッドのテスト"""
        properties = service.calculate_properties(mock_compound)
        
        # プロパティのリストが返されることを確認
        assert isinstance(properties, list)
        assert len(properties) > 0
        
        # 各プロパティがMoleculePropertyインスタンスであることを確認
        for prop in properties:
            assert isinstance(prop, MoleculeProperty)
        
        # サンプル実装では決まったプロパティが返されるはず
        property_names = [p.name for p in properties]
        assert "molecular_weight" in property_names
        assert "logp" in property_names
        assert "h_bond_donors" in property_names
        assert "h_bond_acceptors" in property_names
        assert "rotatable_bonds" in property_names
    
    def test_add_hydrogens_compound(self, service, mock_compound):
        """add_hydrogens()メソッドのテスト（化合物）"""
        # デフォルトpHでの水素原子追加
        compound_with_h = service.add_hydrogens(mock_compound)
        
        # 同じインスタンスが返されることを確認
        assert compound_with_h is mock_compound
        
        # メタデータが追加されていることを確認
        mock_compound.add_metadata.assert_called_once()
        args, _ = mock_compound.add_metadata.call_args
        metadata = args[0]
        assert "hydrogens_added" in metadata
        assert metadata["hydrogens_added"] is True
        assert metadata["ph"] == 7.4  # デフォルト値
        
        # 特定のpH値を指定した場合
        mock_compound.add_metadata.reset_mock()
        
        compound_with_h = service.add_hydrogens(mock_compound, ph=6.5)
        
        mock_compound.add_metadata.assert_called_once()
        args, _ = mock_compound.add_metadata.call_args
        metadata = args[0]
        assert metadata["hydrogens_added"] is True
        assert metadata["ph"] == 6.5  # 指定した値
    
    def test_add_hydrogens_protein(self, service, mock_protein):
        """add_hydrogens()メソッドのテスト（タンパク質）"""
        # デフォルトpHでの水素原子追加
        protein_with_h = service.add_hydrogens(mock_protein)
        
        # 同じインスタンスが返されることを確認
        assert protein_with_h is mock_protein
        
        # メタデータが追加されていることを確認
        mock_protein.add_metadata.assert_called_once()
        args, _ = mock_protein.add_metadata.call_args
        metadata = args[0]
        assert "hydrogens_added" in metadata
        assert metadata["hydrogens_added"] is True
        assert metadata["ph"] == 7.4  # デフォルト値
    
    def test_remove_water(self, service, mock_protein):
        """remove_water()メソッドのテスト"""
        protein_no_water = service.remove_water(mock_protein)
        
        # 同じインスタンスが返されることを確認
        assert protein_no_water is mock_protein
        
        # メタデータが追加されていることを確認
        mock_protein.add_metadata.assert_called_once()
        args, _ = mock_protein.add_metadata.call_args
        metadata = args[0]
        assert "water_removed" in metadata
        assert metadata["water_removed"] is True
    
    def test_get_supported_methods(self, service):
        """get_supported_methods()メソッドのテスト"""
        methods = service.get_supported_methods()
        
        # メソッドの辞書が返されることを確認
        assert isinstance(methods, dict)
        assert len(methods) > 0
        
        # サポートされているメソッドが含まれることを確認
        assert "default" in methods
        assert "minimize" in methods
        assert "keepwater" in methods
        
        # 各メソッドの説明が文字列であることを確認
        for description in methods.values():
            assert isinstance(description, str)
            assert len(description) > 0