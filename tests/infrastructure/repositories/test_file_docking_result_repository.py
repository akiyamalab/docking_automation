"""
FileDockingResultRepositoryのテストモジュール。

このモジュールは、FileDockingResultRepositoryクラスの機能をテストします。
"""

import os
import tempfile
from pathlib import Path
from typing import Generator

import pytest

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.file_docking_result_repository import (
    FileDockingResultRepository,
)


@pytest.fixture
def temp_dir() -> Generator[Path, None, None]:
    """
    一時ディレクトリを作成するフィクスチャ。

    Yields:
        一時ディレクトリのパス
    """
    with tempfile.TemporaryDirectory() as temp_dir:
        yield Path(temp_dir)


@pytest.fixture
def repository(temp_dir: Path) -> FileDockingResultRepository:
    """
    テスト用のFileDockingResultRepositoryを作成するフィクスチャ。

    Args:
        temp_dir: 一時ディレクトリのパス

    Returns:
        テスト用のFileDockingResultRepository
    """
    return FileDockingResultRepository(base_directory=temp_dir)


@pytest.fixture
def sample_result() -> Generator[DockingResult, None, None]:
    """
    サンプルのドッキング結果を作成するフィクスチャ。

    Yields:
        サンプルのドッキング結果
    """
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sdf") as temp_file:
        # サンプルのSDFデータを書き込む
        temp_file.write(
            """
HEADER    PROTEIN: protein1
COMPND    COMPOUND: compound_set1_0
REMARK    DOCKING_SCORE: -8.5
ATOM      1  C   UNK     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   UNK     1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   UNK     1       1.000   1.000   0.000  1.00  0.00           C
ATOM      4  C   UNK     1       0.000   1.000   0.000  1.00  0.00           C
CONECT    1    2    4
CONECT    2    1    3
CONECT    3    2    4
CONECT    4    1    3
END
"""
        )

    # DockingResultオブジェクトを作成
    result = DockingResult.create(
        result_path=Path(temp_file.name),
        protein_id="protein1",
        compound_set_id="compound_set1",
        compound_index=0,
        docking_score=-8.5,
        metadata={
            "compound_name": "Compound_1",
            "binding_mode": "active",
            "interactions": ["hydrogen_bond", "hydrophobic"],
        },
    )

    yield result

    # 一時ファイルを削除
    os.unlink(temp_file.name)


@pytest.fixture
def sample_results() -> Generator[list[DockingResult], None, None]:
    """
    複数のサンプルドッキング結果を作成するフィクスチャ。

    Yields:
        サンプルのドッキング結果のリスト
    """
    temp_files = []
    results = []

    # サンプルデータ
    sample_data = [
        ("protein1", "compound_set1", 0, -8.5),
        ("protein1", "compound_set1", 1, -7.2),
        ("protein1", "compound_set2", 0, -9.1),
        ("protein2", "compound_set1", 0, -6.8),
        ("protein2", "compound_set2", 1, -8.0),
    ]

    for protein_id, compound_set_id, compound_index, docking_score in sample_data:
        # 一時ファイルを作成
        temp_file = tempfile.NamedTemporaryFile(delete=False, mode="w", suffix=".sdf")
        temp_files.append(temp_file.name)

        # サンプルのSDFデータを書き込む
        temp_file.write(
            f"""
HEADER    PROTEIN: {protein_id}
COMPND    COMPOUND: {compound_set_id}_{compound_index}
REMARK    DOCKING_SCORE: {docking_score}
ATOM      1  C   UNK     1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C   UNK     1       1.000   0.000   0.000  1.00  0.00           C
ATOM      3  C   UNK     1       1.000   1.000   0.000  1.00  0.00           C
ATOM      4  C   UNK     1       0.000   1.000   0.000  1.00  0.00           C
CONECT    1    2    4
CONECT    2    1    3
CONECT    3    2    4
CONECT    4    1    3
END
"""
        )
        temp_file.close()

        # DockingResultオブジェクトを作成
        result = DockingResult.create(
            result_path=Path(temp_file.name),
            protein_id=protein_id,
            compound_set_id=compound_set_id,
            compound_index=compound_index,
            docking_score=docking_score,
            metadata={
                "compound_name": f"Compound_{compound_set_id}_{compound_index}",
                "binding_mode": "active",
                "interactions": ["hydrogen_bond", "hydrophobic"],
            },
        )

        results.append(result)

    yield results

    # 一時ファイルを削除
    for temp_file in temp_files:
        os.unlink(temp_file)


class TestFileDockingResultRepository:
    """
    FileDockingResultRepositoryクラスのテスト。
    """

    def test_save_and_load(self, repository: FileDockingResultRepository, sample_result: DockingResult) -> None:
        """
        save()とload()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_result: サンプルのドッキング結果
        """
        # 結果を保存
        repository.save(sample_result)

        # 結果を読み込み
        loaded_result = repository.load(sample_result.id)

        # 読み込んだ結果が正しいことを確認
        assert loaded_result is not None
        assert loaded_result.id == sample_result.id
        assert loaded_result.protein_id == sample_result.protein_id
        assert loaded_result.compound_set_id == sample_result.compound_set_id
        assert loaded_result.compound_index == sample_result.compound_index
        assert loaded_result.docking_score == sample_result.docking_score
        assert loaded_result.metadata == sample_result.metadata
        assert loaded_result.version == sample_result.version + 1  # バージョンがインクリメントされる

    def test_find_by_protein(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        find_by_protein()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # タンパク質IDによる検索
        protein_id = "protein1"
        results = repository.find_by_protein(protein_id)

        # 検索結果が正しいことを確認
        assert len(results) == 3
        for result in results:
            assert result.protein_id == protein_id

    def test_find_by_compound_set(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        find_by_compound_set()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # 化合物セットIDによる検索
        compound_set_id = "compound_set1"
        results = repository.find_by_compound_set(compound_set_id)

        # 検索結果が正しいことを確認
        assert len(results) == 3
        for result in results:
            assert result.compound_set_id == compound_set_id

    def test_find_by_compound(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        find_by_compound()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # 化合物による検索
        compound_set_id = "compound_set1"
        compound_index = 0
        results = repository.find_by_compound(compound_set_id, compound_index)

        # 検索結果が正しいことを確認
        assert len(results) == 2
        for result in results:
            assert result.compound_set_id == compound_set_id
            assert result.compound_index == compound_index

    def test_find_by_protein_and_compound(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        find_by_protein_and_compound()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # タンパク質IDと化合物による検索
        protein_id = "protein1"
        compound_set_id = "compound_set1"
        compound_index = 0
        results = repository.find_by_protein_and_compound(protein_id, compound_set_id, compound_index)

        # 検索結果が正しいことを確認
        assert len(results) == 1
        result = results.get_all()[0]
        assert result.protein_id == protein_id
        assert result.compound_set_id == compound_set_id
        assert result.compound_index == compound_index

    def test_save_collection(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        save_collection()メソッドのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # コレクションを作成
        collection = DockingResultCollection()
        for result in sample_results:
            collection.add(result)

        # コレクションを保存
        repository.save_collection(collection)

        # 保存されたコレクションを検索
        protein_id = "protein1"
        results = repository.find_by_protein(protein_id)

        # 検索結果が正しいことを確認
        assert len(results) == 3
        for result in results:
            assert result.protein_id == protein_id

    def test_transaction(self, repository: FileDockingResultRepository, sample_results: list[DockingResult]) -> None:
        """
        トランザクション管理のテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # トランザクションを開始
        repository.begin_transaction()

        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # トランザクションをコミット
        repository.commit_transaction()

        # 保存されたコレクションを検索
        protein_id = "protein1"
        results = repository.find_by_protein(protein_id)

        # 検索結果が正しいことを確認
        assert len(results) == 3
        for result in results:
            assert result.protein_id == protein_id

    def test_transaction_rollback(
        self, repository: FileDockingResultRepository, sample_results: list[DockingResult]
    ) -> None:
        """
        トランザクションのロールバックのテスト。

        Args:
            repository: テスト用のリポジトリ
            sample_results: サンプルのドッキング結果のリスト
        """
        # トランザクションを開始
        repository.begin_transaction()

        # 結果を保存
        for result in sample_results:
            repository.save(result)

        # トランザクションをロールバック
        repository.rollback_transaction()

        # 保存されたコレクションを検索
        protein_id = "protein1"
        results = repository.find_by_protein(protein_id)

        # 検索結果が空であることを確認（ロールバックされたため）
        assert len(results) == 0
