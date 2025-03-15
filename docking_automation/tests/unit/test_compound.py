#!/usr/bin/env python3
"""
Compoundクラスのテスト

このモジュールは化合物を表すCompoundクラスをテストします。
Compoundクラスはドッキング計算で使用されるリガンドを表すエンティティです。
"""

import os
import pytest
from docking_automation.molecule.entity.compound import Compound


class TestCompound:
    """Compoundクラスのテストケース"""
    
    def test_initialization(self):
        """初期化と基本プロパティのテスト"""
        # 最小限の引数での初期化
        compound = Compound(
            id="C001",
            path="compounds/aspirin.sdf"
        )
        
        assert compound.id == "C001"
        assert compound.path == "compounds/aspirin.sdf"
        assert compound.structure is None
        assert compound.format is None
        assert compound.metadata == {}
        assert compound.name == "aspirin"  # ファイル名から自動設定
        
        # 全オプション付きでの初期化
        structure_mock = {"atoms": [{"element": "C", "x": 0.0, "y": 0.0, "z": 0.0}]}
        format_mock = {"type": "sdf", "version": "1.0"}
        metadata = {"mw": 180.16, "formula": "C9H8O4"}
        
        compound = Compound(
            id="C002",
            path="compounds/ibuprofen.mol2",
            structure=structure_mock,
            format=format_mock,
            metadata=metadata
        )
        
        assert compound.id == "C002"
        assert compound.path == "compounds/ibuprofen.mol2"
        assert compound.structure is structure_mock
        assert compound.format is format_mock
        assert compound.metadata == metadata
        assert compound.name == "ibuprofen"  # ファイル名から自動設定
    
    def test_name_property(self):
        """name属性のゲッター/セッターのテスト"""
        compound = Compound(
            id="C001",
            path="compounds/aspirin.sdf"
        )
        
        # 初期値の確認（ファイル名から自動設定）
        assert compound.name == "aspirin"
        
        # 名前の変更
        compound.name = "acetylsalicylic_acid"
        assert compound.name == "acetylsalicylic_acid"
        
        # 空文字列も許容される
        compound.name = ""
        assert compound.name == ""
    
    def test_add_metadata(self):
        """メタデータ追加のテスト"""
        compound = Compound(
            id="C001",
            path="compounds/aspirin.sdf"
        )
        
        # 初期状態では空のメタデータ
        assert compound.metadata == {}
        
        # メタデータの追加
        compound.add_metadata({"molecular_weight": 180.16})
        assert "molecular_weight" in compound.metadata
        assert compound.metadata["molecular_weight"] == 180.16
        
        # 追加のメタデータ（既存のメタデータは保持される）
        compound.add_metadata({"logP": 1.43, "formula": "C9H8O4"})
        assert "molecular_weight" in compound.metadata  # 既存のメタデータは保持
        assert "logP" in compound.metadata
        assert "formula" in compound.metadata
        assert compound.metadata["logP"] == 1.43
        assert compound.metadata["formula"] == "C9H8O4"
        
        # 既存のキーがある場合は上書き
        compound.add_metadata({"molecular_weight": 180.2})
        assert compound.metadata["molecular_weight"] == 180.2
    
    def test_str_representation(self):
        """文字列表現のテスト"""
        compound = Compound(
            id="C001",
            path="compounds/aspirin.sdf"
        )
        
        # __str__メソッドのテスト
        assert str(compound) == "Compound(C001, name=aspirin)"
        
        # 名前を変更した場合
        compound.name = "acetylsalicylic_acid"
        assert str(compound) == "Compound(C001, name=acetylsalicylic_acid)"
    
    def test_with_various_file_paths(self):
        """様々なファイルパスでのテスト"""
        # 相対パス
        compound1 = Compound(id="C001", path="aspirin.sdf")
        assert compound1.name == "aspirin"
        
        # 絶対パス
        compound2 = Compound(id="C002", path="/path/to/ibuprofen.mol2")
        assert compound2.name == "ibuprofen"
        
        # 拡張子なし
        compound3 = Compound(id="C003", path="compounds/paracetamol")
        assert compound3.name == "paracetamol"
        
        # 複数の拡張子
        compound4 = Compound(id="C004", path="compounds/caffeine.sdf.gz")
        assert compound4.name == "caffeine.sdf"  # 最後の拡張子のみが除去される
        
        # Windowsスタイルのパス
        compound5 = Compound(id="C005", path="C:\\compounds\\naproxen.pdb")
        assert compound5.name == "naproxen"
    
    def test_metadata_none(self):
        """メタデータがNoneの場合のテスト"""
        # メタデータにNoneを指定
        compound = Compound(
            id="C001",
            path="compounds/aspirin.sdf",
            metadata=None
        )
        
        # 空の辞書に初期化されること
        assert compound.metadata == {}
        
        # メタデータの追加が正常に機能すること
        compound.add_metadata({"test_key": "test_value"})
        assert compound.metadata == {"test_key": "test_value"}