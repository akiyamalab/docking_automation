from pathlib import Path

import pytest

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepository,
)


class TestDockingResultRepository:
    """DockingResultRepositoryクラスのテスト"""

    @pytest.fixture
    def repository(self, tmp_path):
        """テスト用のDockingResultRepositoryインスタンスを作成する"""
        # テスト用のリポジトリディレクトリを作成
        repo_dir = tmp_path / "docking_results"
        repo_dir.mkdir()
        return DockingResultRepository()

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

    @pytest.fixture
    def multiple_results(self, tmp_path):
        """複数のテスト用DockingResultインスタンスを作成する"""
        results = []

        # 複数のタンパク質と化合物の組み合わせに対する結果を作成
        for protein_idx in range(2):
            for compound_set_idx in range(2):
                for compound_idx in range(2):
                    result_path = tmp_path / f"result_p{protein_idx}_cs{compound_set_idx}_c{compound_idx}.sdf"
                    result_path.touch()

                    result = DockingResult(
                        result_path=result_path,
                        protein_id=f"protein{protein_idx}",
                        compound_set_id=f"compound_set{compound_set_idx}",
                        compound_index=compound_idx,
                        docking_score=-10.0 + protein_idx + compound_set_idx + compound_idx,
                    )
                    results.append(result)

        return results

    @pytest.mark.skip(reason="未実装のテスト")
    def test_save_and_load(self, repository, sample_result):
        """保存と読み込みのテスト"""
        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_find_by_protein(self, repository, multiple_results):
        """タンパク質IDによる検索のテスト"""
        # 複数の結果を保存
        for result in multiple_results:
            repository.save(result)

        # protein0に関連する結果を検索
        results = repository.find_by_protein("protein0")

        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_find_by_compound_set(self, repository, multiple_results):
        """化合物セットIDによる検索のテスト"""
        # 複数の結果を保存
        for result in multiple_results:
            repository.save(result)

        # compound_set0に関連する結果を検索
        results = repository.find_by_compound_set("compound_set0")

        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_find_by_compound(self, repository, multiple_results):
        """化合物（化合物セットIDと化合物インデックスの組）による検索のテスト"""
        # 複数の結果を保存
        for result in multiple_results:
            repository.save(result)

        # compound_set0の0番目の化合物に関連する結果を検索
        results = repository.find_by_compound("compound_set0", 0)

        pass

    @pytest.mark.skip(reason="未実装のテスト")
    def test_find_by_protein_and_compound(self, repository, multiple_results):
        """タンパク質IDと化合物（化合物セットIDと化合物インデックスの組）による検索のテスト"""
        # 複数の結果を保存
        for result in multiple_results:
            repository.save(result)

        # protein0とcompound_set0の0番目の化合物に関連する結果を検索
        results = repository.find_by_protein_and_compound("protein0", "compound_set0", 0)

        pass
        raise NotImplementedError()
