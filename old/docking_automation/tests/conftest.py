#!/usr/bin/env python3
"""
テスト環境のセットアップを行うconftest.py

このファイルはpytestがテスト実行前に読み込み、テスト環境のセットアップを行います。
"""

import os
import sys
import pytest
from pathlib import Path

# プロジェクトのルートディレクトリをPATHに追加
@pytest.fixture(scope="session", autouse=True)
def setup_path():
    """テスト実行時にプロジェクトルートをパスに追加"""
    # conftest.pyの場所から2階層上がプロジェクトルート
    project_root = Path(__file__).parent.parent.parent.absolute()
    sys.path.insert(0, str(project_root))
    return project_root

# テスト用の一時ディレクトリを提供するフィクスチャ
@pytest.fixture
def temp_dir(tmp_path):
    """テスト用の一時ディレクトリを提供"""
    return tmp_path

# モック化合物データを提供するフィクスチャ
@pytest.fixture
def mock_compound():
    """テスト用のモック化合物データを提供"""
    from docking_automation.molecule.entity.compound import Compound
    
    return Compound(
        id="test_compound",
        path="mock/compound.sdf"
    )

# モックタンパク質データを提供するフィクスチャ
@pytest.fixture
def mock_protein():
    """テスト用のモックタンパク質データを提供"""
    from docking_automation.molecule.entity.protein import Protein
    
    return Protein(
        id="test_protein",
        path="mock/protein.pdb",
        chains={"A"}
    )