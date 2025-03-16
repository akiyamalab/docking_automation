"""
pytestの共通設定ファイル
"""

import pytest
import os
import sys

# プロジェクトのルートディレクトリをPYTHONPATHに追加
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))