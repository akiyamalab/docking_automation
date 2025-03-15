#!/usr/bin/env python3
"""
Proteinクラスのテスト

このモジュールはタンパク質を表すProteinクラスをテストします。
Proteinクラスはドッキング計算で使用されるレセプターを表すエンティティです。
"""

import os
import pytest
from docking_automation.molecule.entity.protein import Protein


class TestProtein:
    """Proteinクラスのテストケース"""
    
    def test_initialization(self):
        """初期化と基本プロパティのテスト"""
        # 最小限の引数での初期化
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        
        assert protein.id == "P001"
        assert protein.path == "proteins/1abc.pdb"
        assert protein.structure is None
        assert protein.format is None
        assert protein.chains == set()  # 空のセット
        assert protein.metadata == {}
        assert protein.name == "1abc"  # ファイル名から自動設定
        
        # チェーン情報付きでの初期化
        chains = {"A", "B", "C"}
        protein = Protein(
            id="P002",
            path="proteins/2xyz.pdb",
            chains=chains
        )
        
        assert protein.id == "P002"
        assert protein.chains == chains
        assert len(protein.chains) == 3
        
        # 全オプション付きでの初期化
        structure_mock = {"atoms": [{"element": "C", "x": 0.0, "y": 0.0, "z": 0.0}]}
        format_mock = {"type": "pdb", "version": "1.0"}
        metadata = {"resolution": 2.5, "temperature": 298.0}
        
        protein = Protein(
            id="P003",
            path="proteins/3def.pdb",
            structure=structure_mock,
            format=format_mock,
            chains={"X", "Y"},
            metadata=metadata
        )
        
        assert protein.id == "P003"
        assert protein.path == "proteins/3def.pdb"
        assert protein.structure is structure_mock
        assert protein.format is format_mock
        assert protein.chains == {"X", "Y"}
        assert protein.metadata == metadata
        assert protein.name == "3def"  # ファイル名から自動設定
    
    def test_name_property(self):
        """name属性のゲッター/セッターのテスト"""
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        
        # 初期値の確認（ファイル名から自動設定）
        assert protein.name == "1abc"
        
        # 名前の変更
        protein.name = "modified_protein"
        assert protein.name == "modified_protein"
        
        # 空文字列も許容される
        protein.name = ""
        assert protein.name == ""
    
    def test_chains(self):
        """チェーン情報のテスト"""
        # chainsが未指定の場合
        protein1 = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        assert protein1.chains == set()  # 空のセット
        
        # 単一チェーンの場合
        protein2 = Protein(
            id="P002",
            path="proteins/2xyz.pdb",
            chains={"A"}
        )
        assert protein2.chains == {"A"}
        assert len(protein2.chains) == 1
        
        # 複数チェーンの場合
        protein3 = Protein(
            id="P003",
            path="proteins/3def.pdb",
            chains={"A", "B", "C", "D"}
        )
        assert protein3.chains == {"A", "B", "C", "D"}
        assert len(protein3.chains) == 4
        
        # 後からチェーン情報を変更
        protein1.chains = {"X", "Y"}
        assert protein1.chains == {"X", "Y"}
        assert len(protein1.chains) == 2
    
    def test_add_metadata(self):
        """メタデータ追加のテスト"""
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        
        # 初期状態では空のメタデータ
        assert protein.metadata == {}
        
        # メタデータの追加
        protein.add_metadata({"resolution": 2.5})
        assert "resolution" in protein.metadata
        assert protein.metadata["resolution"] == 2.5
        
        # 追加のメタデータ（既存のメタデータは保持される）
        protein.add_metadata({"temperature": 298.0, "ph": 7.4})
        assert "resolution" in protein.metadata  # 既存のメタデータは保持
        assert "temperature" in protein.metadata
        assert "ph" in protein.metadata
        assert protein.metadata["temperature"] == 298.0
        assert protein.metadata["ph"] == 7.4
        
        # 既存のキーがある場合は上書き
        protein.add_metadata({"resolution": 2.0})
        assert protein.metadata["resolution"] == 2.0
    
    def test_has_active_site_defined(self):
        """アクティブサイト定義確認のテスト"""
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        
        # 初期状態ではアクティブサイトなし
        assert protein.has_active_site_defined() is False
        
        # アクティブサイト情報を追加
        protein.add_metadata({
            "active_site": {
                "x": 10.0,
                "y": 20.0,
                "z": 30.0,
                "radius": 5.0
            }
        })
        
        # アクティブサイトが定義されていることを確認
        assert protein.has_active_site_defined() is True
        
        # 異なるフォーマットのアクティブサイト情報
        protein2 = Protein(
            id="P002",
            path="proteins/2xyz.pdb",
            metadata={
                "active_site": [
                    {"residue": "A:123", "atom": "CA"},
                    {"residue": "A:124", "atom": "CA"},
                ]
            }
        )
        
        assert protein2.has_active_site_defined() is True
    
    def test_str_representation(self):
        """文字列表現のテスト"""
        # チェーンなし
        protein1 = Protein(
            id="P001",
            path="proteins/1abc.pdb"
        )
        assert str(protein1) == "Protein(P001, name=1abc, chains=0)"
        
        # 単一チェーン
        protein2 = Protein(
            id="P002",
            path="proteins/2xyz.pdb",
            chains={"A"}
        )
        assert str(protein2) == "Protein(P002, name=2xyz, chains=1)"
        
        # 複数チェーン
        protein3 = Protein(
            id="P003",
            path="proteins/3def.pdb",
            chains={"X", "Y", "Z"}
        )
        assert str(protein3) == "Protein(P003, name=3def, chains=3)"
        
        # 名前を変更した場合
        protein3.name = "modified_name"
        assert str(protein3) == "Protein(P003, name=modified_name, chains=3)"
    
    def test_with_various_file_paths(self):
        """様々なファイルパスでのテスト"""
        # 相対パス
        protein1 = Protein(id="P001", path="1abc.pdb")
        assert protein1.name == "1abc"
        
        # 絶対パス
        protein2 = Protein(id="P002", path="/path/to/2xyz.pdb")
        assert protein2.name == "2xyz"
        
        # 拡張子なし
        protein3 = Protein(id="P003", path="proteins/3def")
        assert protein3.name == "3def"
        
        # 複数の拡張子
        protein4 = Protein(id="P004", path="proteins/4ghi.pdb.gz")
        assert protein4.name == "4ghi.pdb"  # 最後の拡張子のみが除去される
        
        # Windowsスタイルのパス
        protein5 = Protein(id="P005", path="C:\\proteins\\5jkl.cif")
        assert protein5.name == "5jkl"
    
    def test_metadata_none(self):
        """メタデータがNoneの場合のテスト"""
        # メタデータにNoneを指定
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb",
            metadata=None
        )
        
        # 空の辞書に初期化されること
        assert protein.metadata == {}
        
        # メタデータの追加が正常に機能すること
        protein.add_metadata({"test_key": "test_value"})
        assert protein.metadata == {"test_key": "test_value"}
    
    def test_chains_none(self):
        """chainsがNoneの場合のテスト"""
        # chainsにNoneを指定
        protein = Protein(
            id="P001",
            path="proteins/1abc.pdb",
            chains=None
        )
        
        # 空のセットに初期化されること
        assert protein.chains == set()
        assert len(protein.chains) == 0