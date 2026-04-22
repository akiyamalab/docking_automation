"""
pytestの共通設定ファイル
"""

import os
import sys

import pytest

# プロジェクトのルートディレクトリをPYTHONPATHに追加
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: mark test as slow-running integration test")
