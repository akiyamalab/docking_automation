from pathlib import Path

import pytest

from docking_automation.molecule.protein import Protein


class TestProtein:
    """Proteinクラスのテスト"""

    @pytest.fixture
    def sample_protein(self, tmp_path):
        """テスト用のProteinインスタンスを作成する"""
        # テスト用のPDBファイルを作成
        pdb_path = tmp_path / "test_protein.pdb"
        pdb_path.write_text("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N")
        return Protein(path=pdb_path, id="test_protein")

    @pytest.mark.skip(reason="未実装のテスト")
    def test_initialization(self, sample_protein):
        """初期化のテスト"""
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_id_assignment(self, tmp_path):
        """IDの割り当てのテスト"""
        # IDを指定した場合
        pdb_path = tmp_path / "test_protein.pdb"
        pdb_path.write_text("ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N")
        protein = Protein(path=pdb_path, id="custom_id")

        # IDを指定しない場合（ファイル名がIDとなる）
        protein_no_id = Protein(path=pdb_path)

        pass
