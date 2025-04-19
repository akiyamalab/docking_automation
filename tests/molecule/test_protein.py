from pathlib import Path

import pytest

from docking_automation.molecule.protein import Protein


class TestProtein:
    """Proteinクラスのテスト"""

    @pytest.fixture
    def create_test_protein(self, tmp_path):
        """テスト用のProteinインスタンスを作成する関数を返すフィクスチャ"""

        def _create(file_name, amino_acid_type="ALA_N", id=None):
            """
            指定されたアミノ酸タイプのProteinインスタンスを作成する

            Args:
                file_name: 作成するファイル名
                amino_acid_type: アミノ酸タイプ
                    - "GLY_N": グリシンのN原子のみ
                    - "ALA_N": アラニンのN原子のみ
                    - "ALA_N_CA": アラニンのN原子とCA原子
                id: Proteinインスタンスに設定するID (Noneの場合はファイル名から自動生成)

            Returns:
                作成されたProteinインスタンス
            """
            # アミノ酸タイプに応じてファイル内容を決定
            if amino_acid_type == "GLY_N":
                content = "ATOM      1  N   GLY A   1      10.000  10.000  10.000  1.00  0.00           N"
            elif amino_acid_type == "ALA_N":
                content = "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N"
            elif amino_acid_type == "ALA_N_CA":
                content = (
                    "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00  0.00           N\n"
                    + "ATOM      2  CA  ALA A   1      11.000  10.000  10.000  1.00  0.00           C"
                )
            else:
                raise ValueError(f"不明なアミノ酸タイプ: {amino_acid_type}")

            # ファイルを作成
            pdb_path = tmp_path / file_name
            pdb_path.write_text(content)

            # Proteinインスタンスを作成して返す
            if id is None:
                # ファイル名から拡張子を除いた部分を取得
                file_stem = Path(file_name).stem
                id = f"test_{file_stem}"
            return Protein(path=pdb_path, id=id)

        return _create

    @pytest.fixture
    def sample_protein(self, create_test_protein):
        """テスト用のProteinインスタンスを作成する"""
        return create_test_protein("test_protein.pdb", amino_acid_type="ALA_N", id="test_protein")

    @pytest.mark.skip(reason="未実装のテスト")
    def test_initialization(self, sample_protein):
        """初期化のテスト"""
        pass

    def test_id_assignment(self, create_test_protein):
        """IDの割り当てのテスト"""
        # IDを指定した場合
        protein = create_test_protein("test_protein.pdb", amino_acid_type="ALA_N", id="custom_id")
        assert protein.id == "custom_id"

        # IDを指定しない場合（ファイル名がIDとなる）
        protein_no_id = create_test_protein("test_protein_no_id.pdb", amino_acid_type="ALA_N")
        assert protein_no_id.id == "test_test_protein_no_id"

    def test_content_hash(self, create_test_protein):
        """ファイル内容に基づいたハッシュ値のテスト"""
        # 同じ内容の2つのProteinインスタンスを作成
        protein1 = create_test_protein("test_protein1.pdb", amino_acid_type="ALA_N", id="test1")
        protein2 = create_test_protein("test_protein2.pdb", amino_acid_type="ALA_N", id="test2")

        # 同じ内容のファイルは同じハッシュ値を持つ
        assert protein1.content_hash == protein2.content_hash

        # ハッシュ値がプロパティに含まれていることを確認
        properties = protein1.get_properties()
        assert "content_hash" in properties
        assert properties["content_hash"] == protein1.content_hash

        # 内容が異なるファイルは異なるハッシュ値を持つ
        protein3 = create_test_protein("test_protein3.pdb", amino_acid_type="GLY_N", id="test3")
        assert protein1.content_hash != protein3.content_hash
