from pathlib import Path

import pytest

from docking_automation.molecule.compound_set import CompoundSet


class TestCompoundSet:
    """CompoundSetクラスのテスト"""

    @pytest.fixture
    def sample_compound_set(self, tmp_path):
        """テスト用のCompoundSetインスタンスを作成する"""
        # テスト用のSDFファイルを作成
        sdf_path = tmp_path / "test_compounds.sdf"
        sdf_content = """
        test_compound1
          RDKit

          9  9  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          1  2  1  0
          2  3  1  0
          3  4  1  0
          4  5  1  0
          5  6  1  0
          6  7  1  0
          7  8  1  0
          8  9  1  0
          9  1  1  0
        M  END
        $$$$
        test_compound2
          RDKit

          9  9  0  0  0  0  0  0  0  0999 V2000
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
            0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0
          1  2  1  0
          2  3  1  0
          3  4  1  0
          4  5  1  0
          5  6  1  0
          6  7  1  0
          7  8  1  0
          8  9  1  0
          9  1  1  0
        M  END
        $$$$
        """
        sdf_path.write_text(sdf_content)
        return CompoundSet(path=sdf_path, id="test_compounds")

    def test_initialization(self, sample_compound_set):
        """初期化のテスト"""
        # 初期化時に指定した値が正しく設定されていることを確認
        assert sample_compound_set.id == "test_compounds"
        assert sample_compound_set.path.name == "test_compounds.sdf"
        assert sample_compound_set.path.exists()

    def test_id_assignment(self, tmp_path):
        """IDの割り当てのテスト"""
        # IDを指定した場合
        sdf_path = tmp_path / "test_compounds.sdf"
        sdf_path.touch()
        compound_set = CompoundSet(path=sdf_path, id="custom_id")
        assert compound_set.id == "custom_id"

        # IDを指定しない場合（ファイル名がIDとなる）
        compound_set_no_id = CompoundSet(path=sdf_path)
        assert compound_set_no_id.id == "test_compounds"

    def test_get_compound_count(self, sample_compound_set):
        """化合物数の取得のテスト"""
        # SDFファイルに含まれる化合物の数を正しく取得できることを確認
        assert sample_compound_set.get_compound_count() == 2

    def test_equality(self, sample_compound_set, tmp_path):
        """等価性比較のテスト"""
        # 同じIDを持つCompoundSetは等価と判断される
        sdf_path = tmp_path / "another_compounds.sdf"
        sdf_path.touch()
        another_compound_set = CompoundSet(path=sdf_path, id="test_compounds")
        assert sample_compound_set == another_compound_set

        # 異なるIDを持つCompoundSetは等価でない
        different_compound_set = CompoundSet(path=sdf_path, id="different_id")
        assert sample_compound_set != different_compound_set

    def test_hash(self, sample_compound_set, tmp_path):
        """ハッシュ値計算のテスト"""
        # 同じIDを持つCompoundSetは同じハッシュ値を持つ
        sdf_path = tmp_path / "another_compounds.sdf"
        sdf_path.touch()
        another_compound_set = CompoundSet(path=sdf_path, id="test_compounds")
        assert hash(sample_compound_set) == hash(another_compound_set)

        # 異なるIDを持つCompoundSetは異なるハッシュ値を持つ
        different_compound_set = CompoundSet(path=sdf_path, id="different_id")
        assert hash(sample_compound_set) != hash(different_compound_set)

    def test_with_new_path(self, sample_compound_set, tmp_path):
        """新しいパスを持つインスタンス作成のテスト"""
        new_path = tmp_path / "new_compounds.sdf"
        new_path.touch()
        new_compound_set = sample_compound_set.with_new_path(new_path)

        # 新しいインスタンスが作成され、パスが更新されていることを確認
        assert new_compound_set.path == new_path
        assert new_compound_set.id == sample_compound_set.id
        assert new_compound_set is not sample_compound_set

        # ドメインイベントが登録されていることを確認
        events = new_compound_set.release_domain_events()
        assert len(events) == 1
        assert events[0].compound_set_id == sample_compound_set.id
        assert events[0].old_path == str(sample_compound_set.path)
        assert events[0].new_path == str(new_path)

    def test_get_compound(self, sample_compound_set):
        """化合物の取得のテスト"""
        # 有効なインデックスの化合物を取得できることを確認
        compound = sample_compound_set.get_compound(0)
        assert compound["index"] == 0
        assert compound["compound_set_id"] == sample_compound_set.id

        # 範囲外のインデックスでIndexErrorが発生することを確認
        with pytest.raises(IndexError):
            sample_compound_set.get_compound(100)

    def test_get_all_compounds(self, sample_compound_set):
        """すべての化合物の取得のテスト"""
        compounds = sample_compound_set.get_all_compounds()
        assert len(compounds) == 2
        assert compounds[0]["index"] == 0
        assert compounds[1]["index"] == 1
        assert compounds[0]["compound_set_id"] == sample_compound_set.id
        assert compounds[1]["compound_set_id"] == sample_compound_set.id

    def test_iteration(self, sample_compound_set):
        """イテレーションのテスト"""
        compounds = list(sample_compound_set)
        assert len(compounds) == 2
        assert compounds[0]["index"] == 0
        assert compounds[1]["index"] == 1

    def test_length(self, sample_compound_set):
        """長さの取得のテスト"""
        assert len(sample_compound_set) == 2

    def test_get_properties(self, sample_compound_set):
        """プロパティの取得のテスト"""
        properties = sample_compound_set.get_properties()
        assert properties["id"] == sample_compound_set.id
        assert properties["path"] == str(sample_compound_set.path)
        assert properties["file_format"] == "sdf"
        assert properties["compound_count"] == 2

    def test_with_index_range(self, sample_compound_set):
        """インデックス範囲を持つインスタンス作成のテスト"""
        # 有効なインデックス範囲でインスタンスを作成
        ranged_compound_set = sample_compound_set.with_index_range(0, 1)
        assert ranged_compound_set.id == f"{sample_compound_set.id}_range_0_1"
        assert len(ranged_compound_set) == 1
        assert ranged_compound_set.get_compound_count() == 2  # 元の化合物数は変わらない

        # 範囲外のインデックスでValueErrorが発生することを確認
        with pytest.raises(ValueError):
            sample_compound_set.with_index_range(-1, 1)
        with pytest.raises(ValueError):
            sample_compound_set.with_index_range(0, 3)
        with pytest.raises(ValueError):
            sample_compound_set.with_index_range(1, 1)

    def test_create_factory_method(self, tmp_path):
        """createファクトリメソッドのテスト"""
        sdf_path = tmp_path / "factory_test.sdf"
        sdf_path.touch()

        # IDを指定した場合
        compound_set = CompoundSet.create(path=sdf_path, id="factory_id")
        assert compound_set.id == "factory_id"
        assert compound_set.path == sdf_path

        # IDを指定しない場合
        compound_set_no_id = CompoundSet.create(path=sdf_path)
        assert compound_set_no_id.id == "factory_test"

        # ドメインイベントが登録されていることを確認
        events = compound_set.release_domain_events()
        assert len(events) == 1
        assert events[0].compound_set_id == compound_set.id
        assert events[0].path == str(sdf_path)

    def test_create_empty_factory_method(self):
        """create_emptyファクトリメソッドのテスト"""
        # 空のCompoundSetを作成
        compound_set = CompoundSet.create_empty(id="empty_test")
        assert compound_set.id == "empty_test"
        assert compound_set.path.exists()
        assert compound_set.path.suffix == ".sdf"
        assert compound_set.get_compound_count() == 0
