import hashlib
import os
from pathlib import Path

import pytest

from docking_automation.infrastructure.utilities.file_utils import (
    calculate_file_content_hash,
)


class TestFileUtils:
    """file_utilsモジュールのテスト"""

    def test_calculate_file_content_hash(self, tmp_path):
        """calculate_file_content_hash関数のテスト"""
        # テスト用のファイルを作成
        test_file_path = tmp_path / "test_file.txt"
        test_content = "This is a test file content."
        test_file_path.write_text(test_content)

        # ハッシュ値を計算
        hash_value = calculate_file_content_hash(test_file_path)

        # 期待されるハッシュ値を計算（直接計算）
        expected_hash = hashlib.sha256(test_content.encode()).hexdigest()

        # ハッシュ値が期待通りであることを確認
        assert hash_value == expected_hash

    def test_calculate_file_content_hash_with_binary_data(self, tmp_path):
        """バイナリデータを含むファイルのハッシュ値計算テスト"""
        # テスト用のバイナリファイルを作成
        test_file_path = tmp_path / "test_binary.bin"
        test_content = bytes([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
        with open(test_file_path, "wb") as f:
            f.write(test_content)

        # ハッシュ値を計算
        hash_value = calculate_file_content_hash(test_file_path)

        # 期待されるハッシュ値を計算（直接計算）
        expected_hash = hashlib.sha256(test_content).hexdigest()

        # ハッシュ値が期待通りであることを確認
        assert hash_value == expected_hash

    def test_calculate_file_content_hash_nonexistent_file(self):
        """存在しないファイルのハッシュ値計算テスト"""
        # 存在しないファイルのパス
        nonexistent_path = Path("nonexistent_file.txt")

        # ハッシュ値を計算（ValueErrorが発生することを期待）
        with pytest.raises(ValueError) as excinfo:
            calculate_file_content_hash(nonexistent_path)

        # エラーメッセージを確認
        assert "指定されたパス" in str(excinfo.value)
        assert "が存在しません" in str(excinfo.value)

    def test_calculate_file_content_hash_large_file(self, tmp_path):
        """大きなファイルのハッシュ値計算テスト"""
        # テスト用の大きなファイルを作成（8KBのファイル）
        test_file_path = tmp_path / "large_file.bin"
        # 8KB（4096バイト * 2）のデータを作成
        test_content = bytes([i % 256 for i in range(4096 * 2)])
        with open(test_file_path, "wb") as f:
            f.write(test_content)

        # ハッシュ値を計算
        hash_value = calculate_file_content_hash(test_file_path)

        # 期待されるハッシュ値を計算（直接計算）
        expected_hash = hashlib.sha256(test_content).hexdigest()

        # ハッシュ値が期待通りであることを確認
        assert hash_value == expected_hash
