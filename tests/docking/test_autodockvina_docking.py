from pathlib import Path

from docking_automation.docking import DockingResult


class TestDockingResult:
    """DockingResultクラスのテスト"""

    def test_docking_result_has_hash_attributes(self):
        """DockingResultオブジェクトがhash属性を持つことを確認する"""
        # テスト用のパラメータ
        result_path = Path("/tmp/test.sdf")
        protein_id = "test_protein"
        compound_set_id = "test_compound_set"
        compound_index = 0
        docking_score = -8.0
        protein_content_hash = "test_protein_hash"
        compoundset_content_hash = "test_compound_hash"
        metadata = {"test": "metadata"}

        # DockingResultオブジェクトを作成
        result = DockingResult(
            result_path=result_path,
            protein_id=protein_id,
            compound_set_id=compound_set_id,
            compound_index=compound_index,
            docking_score=docking_score,
            protein_content_hash=protein_content_hash,
            compoundset_content_hash=compoundset_content_hash,
            metadata=metadata,
        )

        # 検証
        assert hasattr(result, "protein_content_hash")
        assert hasattr(result, "compoundset_content_hash")
        assert result.protein_content_hash == protein_content_hash
        assert result.compoundset_content_hash == compoundset_content_hash
