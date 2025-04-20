from pathlib import Path

import pytest

from docking_automation.docking.docking_result import DockingResult


class TestDockingResult:
    """DockingResultクラスのテスト"""

    @pytest.fixture
    def sample_result(self, tmp_path):
        """テスト用のDockingResultインスタンスを作成する"""
        result_path = tmp_path / "docking_result.sdf"
        result_path.touch()
        return DockingResult(
            result_path=result_path,
            protein_id="protein1",
            compound_set_id="compounds1",
            compound_index=0,
            docking_score=-8.5,
            protein_content_hash="protein_hash_1",
            compoundset_content_hash="compound_hash_1",
            metadata={"tool": "AutoDock Vina"},
        )

    def test_initialization(self, sample_result):
        """初期化のテスト"""
        # 初期化時に指定した値が正しく設定されていることを確認
        assert sample_result.protein_id == "protein1"
        assert sample_result.compound_set_id == "compounds1"
        assert sample_result.compound_index == 0
        assert sample_result.docking_score == -8.5
        assert sample_result.protein_content_hash == "protein_hash_1"
        assert sample_result.compoundset_content_hash == "compound_hash_1"
        assert sample_result.metadata == {"tool": "AutoDock Vina"}
        assert sample_result.id.startswith("protein1_compounds1_0_")
        assert sample_result.version == 1

    def test_get_score(self, sample_result):
        """スコアの取得のテスト"""
        # get_scoreメソッドが正しくドッキングスコアを返すことを確認
        assert sample_result.get_score() == -8.5
        # 直接アクセスと同じ値が返されることを確認
        assert sample_result.get_score() == sample_result.docking_score

    def test_get_pose(self, sample_result):
        """ポーズの取得のテスト"""
        # get_poseメソッドが正しくポーズファイルのパスを返すことを確認
        assert sample_result.get_pose() == sample_result.result_path
        # 返されるパスが存在することを確認
        assert sample_result.get_pose().exists()

    def test_metadata_access(self, sample_result):
        """メタデータへのアクセスのテスト"""
        # get_metadata_valueメソッドが正しくメタデータの値を返すことを確認
        assert sample_result.get_metadata_value("tool") == "AutoDock Vina"
        # 存在しないキーの場合はデフォルト値が返されることを確認
        assert sample_result.get_metadata_value("non_existent_key") is None
        assert sample_result.get_metadata_value("non_existent_key", "default") == "default"

    def test_with_metadata(self, sample_result):
        """メタデータを追加した新しいインスタンス作成のテスト"""
        # 新しいメタデータを持つ新しいインスタンスが作成されることを確認
        new_result = sample_result.with_metadata("new_key", "new_value")

        # 元のインスタンスは変更されていないことを確認
        assert "new_key" not in sample_result.metadata

        # 新しいインスタンスには新しいメタデータが追加されていることを確認
        assert new_result.metadata["new_key"] == "new_value"
        # 元のメタデータも保持されていることを確認
        assert new_result.metadata["tool"] == "AutoDock Vina"

        # 他の属性も同じ値が設定されていることを確認
        assert new_result.protein_id == sample_result.protein_id
        assert new_result.compound_set_id == sample_result.compound_set_id
        assert new_result.compound_index == sample_result.compound_index
        assert new_result.docking_score == sample_result.docking_score
        assert new_result.protein_content_hash == sample_result.protein_content_hash
        assert new_result.compoundset_content_hash == sample_result.compoundset_content_hash
        assert new_result.result_path == sample_result.result_path

    def test_domain_events(self, sample_result):
        """ドメインイベントの登録と取得のテスト"""
        from docking_automation.docking.docking_result_events import (
            DockingResultCreatedEvent,
        )

        # ドメインイベントを登録
        event = DockingResultCreatedEvent(result=sample_result)
        sample_result.register_domain_event(event)

        # 登録されたイベントを取得
        events = sample_result.get_domain_events()
        assert len(events) == 1
        assert event in events

        # イベントをクリア
        sample_result.clear_domain_events()
        assert len(sample_result.get_domain_events()) == 0

    def test_with_domain_event(self, sample_result):
        """ドメインイベントを含む新しいインスタンス作成のテスト"""
        from docking_automation.docking.docking_result_events import (
            DockingResultCreatedEvent,
        )

        # ドメインイベントを含む新しいインスタンスを作成
        event = DockingResultCreatedEvent(result=sample_result)
        new_result = sample_result.with_domain_event(event)

        # 元のインスタンスは変更されていないことを確認
        assert len(sample_result.get_domain_events()) == 0

        # 新しいインスタンスにはイベントが登録されていることを確認
        events = new_result.get_domain_events()
        assert len(events) == 1
        assert event in events

    def test_with_incremented_version(self, sample_result):
        """バージョンをインクリメントした新しいインスタンス作成のテスト"""
        # バージョンをインクリメントした新しいインスタンスを作成
        new_result = sample_result.with_incremented_version()

        # 元のインスタンスは変更されていないことを確認
        assert sample_result.version == 1

        # 新しいインスタンスのバージョンがインクリメントされていることを確認
        assert new_result.version == 2

        # 他の属性は同じ値が設定されていることを確認
        assert new_result.protein_id == sample_result.protein_id
        assert new_result.compound_set_id == sample_result.compound_set_id
        assert new_result.compound_index == sample_result.compound_index
        assert new_result.docking_score == sample_result.docking_score
        assert new_result.protein_content_hash == sample_result.protein_content_hash
        assert new_result.compoundset_content_hash == sample_result.compoundset_content_hash
        assert new_result.result_path == sample_result.result_path
        assert new_result.metadata == sample_result.metadata

    def test_equality(self, sample_result, tmp_path):
        """等価性比較のテスト"""
        # 同じIDを持つDockingResultは等価と判断される
        result_path = tmp_path / "another_result.sdf"
        result_path.touch()
        another_result = DockingResult(
            result_path=result_path,
            protein_id="protein2",
            compound_set_id="compounds2",
            compound_index=1,
            docking_score=-7.5,
            protein_content_hash="protein_hash_2",
            compoundset_content_hash="compound_hash_2",
            id=sample_result.id,
        )
        assert sample_result == another_result

        # 異なるIDを持つDockingResultは等価でない
        different_result = DockingResult(
            result_path=result_path,
            protein_id="protein2",
            compound_set_id="compounds2",
            compound_index=1,
            docking_score=-7.5,
            protein_content_hash="protein_hash_2",
            compoundset_content_hash="compound_hash_2",
        )
        assert sample_result != different_result

    def test_hash(self, sample_result, tmp_path):
        """ハッシュ値計算のテスト"""
        # 同じIDを持つDockingResultは同じハッシュ値を持つ
        result_path = tmp_path / "another_result.sdf"
        result_path.touch()
        another_result = DockingResult(
            result_path=result_path,
            protein_id="protein2",
            compound_set_id="compounds2",
            compound_index=1,
            docking_score=-7.5,
            protein_content_hash="protein_hash_2",
            compoundset_content_hash="compound_hash_2",
            id=sample_result.id,
        )
        assert hash(sample_result) == hash(another_result)

        # 異なるIDを持つDockingResultは異なるハッシュ値を持つ
        different_result = DockingResult(
            result_path=result_path,
            protein_id="protein2",
            compound_set_id="compounds2",
            compound_index=1,
            docking_score=-7.5,
            protein_content_hash="protein_hash_2",
            compoundset_content_hash="compound_hash_2",
        )
        assert hash(sample_result) != hash(different_result)

    def test_format_result(self, sample_result):
        """結果の整形のテスト"""
        # ランクなしの場合
        formatted = sample_result.format_result()
        assert "Score: -8.5" in formatted
        assert "Compound: actives_subset_0" in formatted

        # ランクありの場合
        formatted = sample_result.format_result(rank=1)
        assert "1. Score: -8.5" in formatted

        # 元のインデックスを指定した場合
        formatted = sample_result.format_result(original_index=10)
        assert "Compound: actives_subset_10" in formatted

    def test_get_compound_info(self, sample_result):
        """化合物情報の取得のテスト"""
        # メタデータに化合物名がない場合
        info = sample_result.get_compound_info()
        assert "化合物ID: actives_subset_0" in info

        # メタデータに化合物名がある場合
        new_result = sample_result.with_metadata("compound_name", "Test Compound")
        info = new_result.get_compound_info()
        assert "化合物名: Test Compound" in info

        # 元のインデックスを指定した場合
        info = sample_result.get_compound_info(original_index=10)
        assert "化合物ID: actives_subset_10" in info

    def test_get_original_index(self, sample_result):
        """元のインデックスの計算のテスト"""
        # 化合物セットのリストが空の場合
        original_index = sample_result.get_original_index([])
        assert original_index == sample_result.compound_index

        # 一致する化合物セットがない場合
        import tempfile

        from docking_automation.molecule.compound_set import CompoundSet

        with tempfile.TemporaryDirectory() as temp_dir:
            sdf_path = Path(temp_dir) / "test.sdf"
            sdf_path.touch()
            compound_set = CompoundSet(path=sdf_path, id="different_id")
            original_index = sample_result.get_original_index([compound_set])
            assert original_index == sample_result.compound_index

        # 一致する化合物セットがあり、インデックス範囲が設定されている場合
        class MockCompoundSet:
            @property
            def id(self):
                return "compounds1"

            def get_properties(self):
                return {"index_range": {"start": 100, "end": 200}}

        mock_compound_set = MockCompoundSet()
        original_index = sample_result.get_original_index([mock_compound_set])
        assert original_index == 100 + sample_result.compound_index

    def test_create_factory_method(self, tmp_path):
        """createファクトリメソッドのテスト"""
        result_path = tmp_path / "factory_result.sdf"
        result_path.touch()

        # ファクトリメソッドでインスタンスを作成
        result = DockingResult.create(
            result_path=result_path,
            protein_id="protein_factory",
            compound_set_id="compounds_factory",
            compound_index=5,
            docking_score=-9.0,
            protein_content_hash="protein_hash_factory",
            compoundset_content_hash="compound_hash_factory",
            metadata={"factory": True},
            id="custom_id",
        )

        # 指定した値が正しく設定されていることを確認
        assert result.protein_id == "protein_factory"
        assert result.compound_set_id == "compounds_factory"
        assert result.compound_index == 5
        assert result.docking_score == -9.0
        assert result.protein_content_hash == "protein_hash_factory"
        assert result.compoundset_content_hash == "compound_hash_factory"
        assert result.metadata == {"factory": True}
        assert result.id == "custom_id"
        assert result.version == 1

        # ドメインイベントが登録されていることを確認
        events = result.get_domain_events()
        assert len(events) == 1
        from docking_automation.docking.docking_result_events import (
            DockingResultCreatedEvent,
        )

        assert isinstance(events.pop(), DockingResultCreatedEvent)
