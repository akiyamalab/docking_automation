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

    @pytest.mark.skip(reason="未実装のテスト")
    def test_initialization(self, sample_result):
        """初期化のテスト"""
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_get_score(self, sample_result):
        """スコアの取得のテスト"""
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_get_pose(self, sample_result):
        """ポーズの取得のテスト"""
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_metadata_access(self, sample_result):
        """メタデータへのアクセスのテスト"""
        pass
