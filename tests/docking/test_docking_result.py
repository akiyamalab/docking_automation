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
            metadata={"tool": "AutoDock Vina"},
        )

    def test_initialization(self, sample_result):
        """初期化のテスト"""
        # 初期化時に指定した値が正しく設定されていることを確認
        assert sample_result.protein_id == "protein1"
        assert sample_result.compound_set_id == "compounds1"
        assert sample_result.compound_index == 0
        assert sample_result.docking_score == -8.5
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
