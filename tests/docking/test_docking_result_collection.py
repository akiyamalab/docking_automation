from pathlib import Path

import pytest

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection


class TestDockingResultCollection:
    """DockingResultCollectionクラスのテスト"""

    @pytest.fixture
    def sample_results(self, tmp_path):
        """テスト用のDockingResultインスタンスのリストを作成する"""
        results = []

        # 複数のタンパク質と化合物の組み合わせに対する結果を作成
        for protein_idx in range(2):
            for compound_set_idx in range(2):
                for compound_idx in range(2):
                    result_path = tmp_path / f"result_p{protein_idx}_cs{compound_set_idx}_c{compound_idx}.sdf"
                    result_path.touch()

                    # スコアは低いほど良い（負の値が大きいほど良い）とする
                    score = -10.0 + protein_idx + compound_set_idx + compound_idx

                    result = DockingResult(
                        result_path=result_path,
                        protein_id=f"protein{protein_idx}",
                        compound_set_id=f"compound_set{compound_set_idx}",
                        compound_index=compound_idx,
                        docking_score=score,
                        protein_content_hash=f"protein_hash_{protein_idx}",
                        compoundset_content_hash=f"compound_hash_{compound_set_idx}",
                        metadata={"score": score},
                    )
                    results.append(result)

        return results

    @pytest.fixture
    def sample_collection(self, sample_results):
        """テスト用のDockingResultCollectionインスタンスを作成する"""
        collection = DockingResultCollection()
        for result in sample_results:
            collection.add(result)
        return collection

    def test_initialization(self):
        """初期化のテスト"""
        # 空のコレクションの初期化
        collection = DockingResultCollection()
        assert len(collection) == 0

        # 結果リストを指定した初期化
        results = [
            DockingResult(
                result_path=Path("/tmp/result1.sdf"),
                protein_id="protein1",
                compound_set_id="compounds1",
                compound_index=0,
                docking_score=-8.5,
                protein_content_hash="protein_hash_test",
                compoundset_content_hash="compound_hash_test",
            ),
            DockingResult(
                result_path=Path("/tmp/result2.sdf"),
                protein_id="protein1",
                compound_set_id="compounds1",
                compound_index=1,
                docking_score=-7.5,
                protein_content_hash="protein_hash_test",
                compoundset_content_hash="compound_hash_test",
            ),
        ]
        collection_with_results = DockingResultCollection(results)
        assert len(collection_with_results) == 2

    def test_add(self, sample_results):
        """結果の追加のテスト"""
        collection = DockingResultCollection()
        # 結果を1つずつ追加
        for result in sample_results:
            collection.add(result)

        # 追加した結果の数が正しいことを確認
        assert len(collection) == len(sample_results)

        # 追加した結果がすべて含まれていることを確認
        for result in sample_results:
            assert result in collection.get_all()

    def test_extend(self, sample_results):
        """複数の結果を一度に追加するテスト"""
        collection = DockingResultCollection()
        half_size = len(sample_results) // 2

        # リストとして結果を追加
        collection.extend(sample_results[:half_size])
        assert len(collection) == half_size

        # コレクションとして結果を追加
        other_collection = DockingResultCollection(sample_results[half_size:])
        collection.extend(other_collection)

        # 追加した結果の数が正しいことを確認
        assert len(collection) == len(sample_results)

        # 追加した結果がすべて含まれていることを確認
        for result in sample_results:
            assert result in collection.get_all()

    def test_get_all(self, sample_collection, sample_results):
        """すべての結果の取得のテスト"""
        # すべての結果を取得
        all_results = sample_collection.get_all()

        # 結果の数が正しいことを確認
        assert len(all_results) == len(sample_results)

        # 結果がスコア順にソートされていることを確認
        for i in range(1, len(all_results)):
            assert all_results[i - 1].docking_score <= all_results[i].docking_score

    def test_get_top(self, sample_collection):
        """上位n件の結果の取得のテスト"""
        # 上位3件の結果を取得
        top_3 = sample_collection.get_top(3)

        # 結果の数が正しいことを確認
        assert len(top_3) == 3

        # 結果がスコア順にソートされていることを確認
        for i in range(1, len(top_3)):
            assert top_3[i - 1].docking_score <= top_3[i].docking_score

        # コレクションサイズより大きい数を指定した場合
        all_results = sample_collection.get_all()
        top_all = sample_collection.get_top(len(sample_collection) + 10)
        assert len(top_all) == len(all_results)

    def test_filter(self, sample_collection):
        """フィルタリングのテスト"""
        # スコアが-9.5より良い（より負の値が大きい）結果のみをフィルタリング
        filtered = sample_collection.filter(lambda result: result.docking_score < -9.5)

        # フィルタリング結果の検証
        for result in filtered:
            assert result.docking_score < -9.5

        # 特定のタンパク質IDに関連する結果のみをフィルタリング
        protein_filtered = sample_collection.filter(lambda result: result.protein_id == "protein0")

        # フィルタリング結果の検証
        for result in protein_filtered:
            assert result.protein_id == "protein0"

        # 特定の化合物に関連する結果のみをフィルタリング
        compound_filtered = sample_collection.filter(
            lambda result: result.compound_set_id == "compound_set0" and result.compound_index == 0
        )

        # フィルタリング結果の検証
        for result in compound_filtered:
            assert result.compound_set_id == "compound_set0"
            assert result.compound_index == 0

    def test_merge(self, sample_results):
        """マージのテスト"""
        half_size = len(sample_results) // 2
        collection1 = DockingResultCollection(sample_results[:half_size])
        collection2 = DockingResultCollection(sample_results[half_size:])

        # コレクションをマージ
        merged = collection1.merge(collection2)

        # マージ結果の検証
        assert len(merged) == len(sample_results)

        # 元のコレクションが変更されていないことを確認
        assert len(collection1) == half_size
        assert len(collection2) == len(sample_results) - half_size

        # マージされたコレクションにすべての結果が含まれていることを確認
        for result in sample_results:
            assert result in merged.get_all()

    def test_sorted_order(self, sample_results):
        """ソート順序のテスト"""
        # 結果をランダムな順序で追加しても、常にスコア順にソートされることを確認
        import random

        shuffled_results = sample_results.copy()
        random.shuffle(shuffled_results)

        collection = DockingResultCollection()
        for result in shuffled_results:
            collection.add(result)

        # 結果がスコア順にソートされていることを確認
        all_results = collection.get_all()
        for i in range(1, len(all_results)):
            assert all_results[i - 1].docking_score <= all_results[i].docking_score
