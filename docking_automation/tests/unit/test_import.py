#!/usr/bin/env python3
"""
docking_automationパッケージのインポートテスト

このテストはdocking_automationパッケージが正しくインポートできることを確認します。
"""

import sys
import os
import pytest
import subprocess

def test_import_package():
    """docking_automationパッケージのインポートをテスト"""
    # パッケージがインポートできることを確認
    import docking_automation
    
    # インポートが成功したことを確認
    assert docking_automation is not None
    assert hasattr(docking_automation, '__file__')
    
    # パスが正しいことを確認
    package_path = os.path.normpath(docking_automation.__file__)
    expected_suffix = os.path.normpath('docking_automation/__init__.py')
    
    assert package_path.endswith(expected_suffix), \
        f"Package path {package_path} does not end with {expected_suffix}"

def test_directory_structure():
    """ディレクトリ構造のチェック"""
    # docking_automationパッケージのルートディレクトリを取得
    import docking_automation
    package_dir = os.path.dirname(docking_automation.__file__)
    
    # 重要なサブディレクトリが存在することを確認
    expected_dirs = ['molecule', 'docking', 'application', 'infrastructure']
    
    for expected_dir in expected_dirs:
        dir_path = os.path.join(package_dir, expected_dir)
        assert os.path.isdir(dir_path), f"Expected directory {expected_dir} not found"
        
        # 各ディレクトリに__init__.pyファイルが存在することを確認
        init_file = os.path.join(dir_path, '__init__.py')
        assert os.path.isfile(init_file), f"__init__.py not found in {expected_dir}"