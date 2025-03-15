#!/usr/bin/env python3
"""
Receptorクラスのテスト

このモジュールはドッキング計算用のレセプター（準備されたタンパク質）を表すReceptorクラスをテストします。
"""

import pytest
from unittest.mock import MagicMock, patch
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.molecule.entity.protein import Protein
from docking_automation.docking.value_object.grid_box import GridBox


class TestReceptor:
    """Receptorクラスのテストケース"""
    
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
    
    @pytest.fixture
    def mock_grid_box(self):
        """テスト用のモックグリッドボックス"""
        return GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
    
    def test_initialization(self, mock_protein, mock_grid_box):
        """初期化と基本プロパティのテスト"""
        # 最小限の引数での初期化（アクティブサイトは自動生成）
        receptor = Receptor(protein=mock_protein)
        
        assert receptor.protein is mock_protein
        assert receptor.metadata == {}
        assert receptor.name == mock_protein.id  # proteinのIDがデフォルトのname
        assert receptor.active_site is not None  # デフォルトのアクティブサイトが生成されている
        
        # アクティブサイト、名前、メタデータを指定した初期化
        metadata = {"resolution": 2.5}
        receptor = Receptor(
            protein=mock_protein,
            active_site=mock_grid_box,
            metadata=metadata,
            name="receptor_1"
        )
        
        assert receptor.protein is mock_protein
        assert receptor.active_site is mock_grid_box
        assert receptor.metadata is metadata
        assert receptor.name == "receptor_1"
    
    def test_name_initialization_from_metadata(self, mock_protein):
        """メタデータからの名前初期化テスト"""
        # タンパク質メタデータに名前がある場合
        mock_protein.metadata = {"name": "protein_in_metadata"}
        
        receptor = Receptor(protein=mock_protein)
        
        # メタデータの名前が優先される
        assert receptor.name == "protein_in_metadata"
        
        # 明示的に名前を指定した場合は、それが優先される
        receptor = Receptor(protein=mock_protein, name="explicit_name")
        
        assert receptor.name == "explicit_name"
    
    def test_is_prepared(self, mock_protein):
        """準備状態確認のテスト"""
        receptor = Receptor(protein=mock_protein)
        
        # Receptorクラスのis_preparedは常にTrueを返す（テスト用簡略化実装）
        assert receptor.is_prepared() is True
    
    def test_get_path(self, mock_protein):
        """パス取得のテスト"""
        receptor = Receptor(protein=mock_protein)
        
        # proteinのパスが返される
        mock_protein.path = "test/path.pdb"
        assert receptor.get_path() == "test/path.pdb"
        
        # proteinのパスがNoneの場合
        mock_protein.path = None
        assert receptor.get_path() is None
    
    def test_get_id(self, mock_protein):
        """ID取得のテスト"""
        receptor = Receptor(protein=mock_protein)
        
        # proteinのIDが返される
        mock_protein.id = "test_id_123"
        assert receptor.get_id() == "test_id_123"
    
    def test_get_chains(self, mock_protein):
        """チェーン取得のテスト"""
        receptor = Receptor(protein=mock_protein)
        
        # proteinのchainsが返される
        mock_protein.chains = {"A", "B", "C"}
        chains = receptor.get_chains()
        
        assert chains == {"A", "B", "C"}
    
    def test_get_suggested_grid_box(self, mock_protein, mock_grid_box):
        """推奨グリッドボックス取得のテスト"""
        # アクティブサイトなしの場合
        receptor_no_active_site = Receptor(
            protein=mock_protein,
            active_site=None
        )
        
        # タンパク質の中心に標準サイズのグリッドボックスが返される
        grid_box = receptor_no_active_site.get_suggested_grid_box()
        
        assert grid_box is not None
        assert grid_box.center_x == 0.0  # テスト用デフォルト値
        assert grid_box.center_y == 0.0
        assert grid_box.center_z == 0.0
        assert grid_box.size_x == 20.0
        assert grid_box.size_y == 20.0
        assert grid_box.size_z == 20.0
        
        # アクティブサイトありの場合
        receptor_with_active_site = Receptor(
            protein=mock_protein,
            active_site=mock_grid_box
        )
        
        # 指定されたアクティブサイトが返される
        grid_box = receptor_with_active_site.get_suggested_grid_box()
        
        assert grid_box is mock_grid_box
    
    def test_has_active_site_defined(self, mock_protein, mock_grid_box):
        """アクティブサイト定義確認のテスト"""
        # アクティブサイトなしの場合
        receptor = Receptor(
            protein=mock_protein,
            active_site=None  # __post_init__でデフォルトのアクティブサイトが設定される
        )
        
        # __post_init__でデフォルトのグリッドボックスが設定されるため、has_active_site_definedはTrue
        assert receptor.has_active_site_defined() is True
        
        # アクティブサイトありの場合
        receptor = Receptor(
            protein=mock_protein,
            active_site=mock_grid_box
        )
        
        assert receptor.has_active_site_defined() is True
    
    def test_with_active_site(self, mock_protein, mock_grid_box):
        """アクティブサイト変更のテスト"""
        # 元のレセプター
        original_receptor = Receptor(
            protein=mock_protein,
            name="original_receptor",
            metadata={"key": "value"}
        )
        
        # 新しいグリッドボックス
        new_grid_box = GridBox(
            center_x=1.0,
            center_y=2.0,
            center_z=3.0,
            size_x=10.0,
            size_y=10.0,
            size_z=10.0
        )
        
        # アクティブサイトを変更した新しいレセプターを作成
        new_receptor = original_receptor.with_active_site(new_grid_box)
        
        # 元のレセプターと異なるインスタンスであることを確認
        assert new_receptor is not original_receptor
        
        # タンパク質、名前、メタデータは同じであることを確認
        assert new_receptor.protein is original_receptor.protein
        assert new_receptor.name == original_receptor.name
        assert new_receptor.metadata == original_receptor.metadata
        
        # アクティブサイトが更新されていることを確認
        assert new_receptor.active_site is new_grid_box
        assert new_receptor.active_site is not original_receptor.active_site
    
    def test_str_representation(self, mock_protein, mock_grid_box):
        """文字列表現のテスト"""
        # アクティブサイトなしの場合
        receptor = Receptor(
            protein=mock_protein,
            name="test_receptor",
            active_site=None  # __post_init__でデフォルトのアクティブサイトが設定される
        )
        
        str_rep = str(receptor)
        assert "test_receptor" in str_rep
        assert "prepared=True" in str_rep
        assert "active_site_defined=True" in str_rep
        
        # アクティブサイトありの場合
        receptor = Receptor(
            protein=mock_protein,
            name="test_receptor",
            active_site=mock_grid_box
        )
        
        str_rep = str(receptor)
        assert "test_receptor" in str_rep
        assert "prepared=True" in str_rep
        assert "active_site_defined=True" in str_rep