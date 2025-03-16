#!/usr/bin/env python3
"""
MoleculePreparationFactoryクラスのテスト

このモジュールは分子準備サービスのファクトリクラスをテストします。
このファクトリは各ドッキングツールに適した分子準備サービスのインスタンスを生成します。
"""

import pytest
from docking_automation.molecule.service.molecule_preparation_factory import MoleculePreparationFactory
from docking_automation.molecule.service.molecule_preparation_service import MoleculePreparationService
from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService


class TestMoleculePreparationFactory:
    """MoleculePreparationFactoryクラスのテストケース"""
    
    def test_create_for_tool_vina(self):
        """Vina用の分子準備サービス生成テスト"""
        # Vina用サービスの生成
        service = MoleculePreparationFactory.create_for_tool("vina")
        
        # 正しいクラスのインスタンスが生成されることを確認
        assert isinstance(service, VinaPreparationService)
        assert isinstance(service, MoleculePreparationService)  # 基底クラスをチェック
    
    def test_create_for_tool_vina_case_insensitive(self):
        """ツール名の大文字小文字を区別しないテスト"""
        # 大文字でのツール名指定
        service1 = MoleculePreparationFactory.create_for_tool("VINA")
        assert isinstance(service1, VinaPreparationService)
        
        # 混合文字でのツール名指定
        service2 = MoleculePreparationFactory.create_for_tool("Vina")
        assert isinstance(service2, VinaPreparationService)
    
    def test_create_for_unsupported_tool(self):
        """未サポートのツール名指定時のエラーテスト"""
        # 存在しないツール名を指定して例外が発生することを確認
        with pytest.raises(ValueError) as excinfo:
            MoleculePreparationFactory.create_for_tool("unsupported_tool")
        
        # エラーメッセージに適切な情報が含まれることを確認
        error_msg = str(excinfo.value)
        assert "未サポートのドッキングツール" in error_msg
        assert "unsupported_tool" in error_msg
        assert "サポートされているツール" in error_msg
        assert "vina" in error_msg  # サポートされているツール名が含まれる
    
    def test_get_supported_tools(self):
        """サポートされているツール一覧取得テスト"""
        # サポートされているツール一覧を取得
        supported_tools = MoleculePreparationFactory.get_supported_tools()
        
        # 戻り値が辞書であり、キーと値が適切であることを確認
        assert isinstance(supported_tools, dict)
        assert "vina" in supported_tools
        assert "AutoDock Vina" in supported_tools["vina"]