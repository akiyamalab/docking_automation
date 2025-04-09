import time
from multiprocessing import Process
from pathlib import Path

import h5py
import pytest

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
    HDF5DockingResultRepository,
)


@pytest.fixture
def hdf5_repo(tmp_path: Path) -> HDF5DockingResultRepository:
    """一時的なHDF5ファイルを使用するリポジトリのfixture。"""
    hdf5_file = tmp_path / "test_results.hdf5"
    return HDF5DockingResultRepository(hdf5_file_path=hdf5_file)


@pytest.fixture
def sample_result1(tmp_path: Path) -> DockingResult:
    """テスト用のサンプルDockingResult 1。"""
    # result_pathは一時ファイルとして作成
    result_file = tmp_path / "proteinA_set1_0.sdf"
    result_file.touch()
    return DockingResult.create(  # ファクトリメソッドを使用
        result_path=result_file,
        protein_id="proteinA",
        compound_set_id="set1",
        compound_index=0,
        docking_score=-9.5,
        metadata={"pose_data": "pose data for compound1"},  # pose_dataをmetadataに含める
    )  # 閉じ括弧を追加


@pytest.fixture
def sample_result2(tmp_path: Path) -> DockingResult:
    """テスト用のサンプルDockingResult 2。"""
    result_file = tmp_path / "proteinA_set1_1.sdf"
    result_file.touch()
    return DockingResult.create(
        result_path=result_file,
        protein_id="proteinA",
        compound_set_id="set1",
        compound_index=1,  # 異なるインデックス
        docking_score=-8.0,
        metadata={"pose_data": "pose data for compound2"},
    )  # 閉じ括弧を追加


@pytest.fixture
def sample_result3(tmp_path: Path) -> DockingResult:
    """テスト用のサンプルDockingResult 3 (異なるタンパク質)。"""
    result_file = tmp_path / "proteinB_set2_0.sdf"
    result_file.touch()
    return DockingResult.create(
        result_path=result_file,
        protein_id="proteinB",  # 異なるタンパク質ID
        compound_set_id="set2",  # 異なるセットID
        compound_index=0,
        docking_score=-10.0,
        metadata={"pose_data": "pose data for compound3"},
    )


# クラス定義のインデントを修正
class TestHDF5DockingResultRepository:
    """HDF5DockingResultRepositoryのテストクラス。"""

    def test_init_creates_directory(self, tmp_path: Path):
        """リポジトリ初期化時にディレクトリが作成されることを確認する。"""
        repo_dir = tmp_path / "repo_test"
        hdf5_file = repo_dir / "results.hdf5"
        assert not repo_dir.exists()
        HDF5DockingResultRepository(hdf5_file_path=hdf5_file)
        assert repo_dir.exists()
        assert repo_dir.is_dir()

    def test_save_and_load(self, hdf5_repo: HDF5DockingResultRepository, sample_result1: DockingResult):
        """単一の結果を保存し、正しくロードできることを確認する。"""
        hdf5_repo.save(sample_result1)

        # loadメソッドの引数を修正
        loaded_result = hdf5_repo.load(
            sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index
        )

        assert loaded_result is not None
        assert loaded_result.protein_id == sample_result1.protein_id
        assert loaded_result.compound_set_id == sample_result1.compound_set_id
        assert loaded_result.compound_index == sample_result1.compound_index
        assert loaded_result.docking_score == sample_result1.docking_score
        assert loaded_result.get_metadata_value("pose_data") == sample_result1.get_metadata_value("pose_data")
        assert loaded_result.result_path == sample_result1.result_path  # result_pathも確認

    def test_load_non_existent(self, hdf5_repo: HDF5DockingResultRepository):
        """存在しない結果をロードしようとするとNoneが返ることを確認する。"""
        # loadメソッドの引数を修正
        loaded_result = hdf5_repo.load("non_existent_protein", "non_existent_set", 0)
        assert loaded_result is None

    def test_save_overwrites_existing(
        self, hdf5_repo: HDF5DockingResultRepository, sample_result1: DockingResult, tmp_path: Path
    ):
        """同じキーで保存すると既存の結果が上書きされることを確認する。"""
        hdf5_repo.save(sample_result1)

        # 更新用のDockingResultを作成
        updated_result_path = tmp_path / "updated_proteinA_set1_0.sdf"
        updated_result_path.touch()
        updated_result = DockingResult.create(
            result_path=updated_result_path,  # パスも更新される可能性がある
            protein_id=sample_result1.protein_id,
            compound_set_id=sample_result1.compound_set_id,
            compound_index=sample_result1.compound_index,
            docking_score=-10.0,  # スコア更新
            metadata={"pose_data": "updated pose data"},  # メタデータ更新
        )
        hdf5_repo.save(updated_result)

        # loadメソッドの引数を修正
        loaded_result = hdf5_repo.load(
            sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index
        )
        assert loaded_result is not None
        assert loaded_result.docking_score == -10.0
        assert loaded_result.get_metadata_value("pose_data") == "updated pose data"
        assert loaded_result.result_path == updated_result_path  # result_pathが更新されていることを確認

    def test_load_all(
        self,
        hdf5_repo: HDF5DockingResultRepository,
        sample_result1: DockingResult,
        sample_result2: DockingResult,
        sample_result3: DockingResult,
    ):
        """複数の結果を保存し、load_allで全てロードできることを確認する。"""
        hdf5_repo.save(sample_result1)
        hdf5_repo.save(sample_result2)
        hdf5_repo.save(sample_result3)

        collection = hdf5_repo.load_all()
        assert isinstance(collection, DockingResultCollection)
        # 内部属性 _results を使用
        assert len(collection._results) == 3

        # 結果の内容を確認 (順序は不定なので、キーで検索)
        # キーを (protein_id, compound_set_id, compound_index) に変更
        results_dict = {(r.protein_id, r.compound_set_id, r.compound_index): r for r in collection._results}
        key1 = (sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index)
        key2 = (sample_result2.protein_id, sample_result2.compound_set_id, sample_result2.compound_index)
        key3 = (sample_result3.protein_id, sample_result3.compound_set_id, sample_result3.compound_index)

        assert key1 in results_dict
        assert results_dict[key1].docking_score == sample_result1.docking_score
        assert key2 in results_dict
        assert results_dict[key2].docking_score == sample_result2.docking_score
        assert key3 in results_dict
        assert results_dict[key3].docking_score == sample_result3.docking_score

    def test_load_all_empty(self, hdf5_repo: HDF5DockingResultRepository):
        """空のリポジトリからload_allを実行すると空のコレクションが返ることを確認する。"""
        collection = hdf5_repo.load_all()
        assert isinstance(collection, DockingResultCollection)
        # 内部属性 _results を使用
        assert len(collection._results) == 0

    def test_delete(
        self, hdf5_repo: HDF5DockingResultRepository, sample_result1: DockingResult, sample_result2: DockingResult
    ):
        """結果を削除できることを確認する。"""
        hdf5_repo.save(sample_result1)
        hdf5_repo.save(sample_result2)

        # load/deleteメソッドの引数を修正
        assert (
            hdf5_repo.load(sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index)
            is not None
        )
        hdf5_repo.delete(sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index)
        assert (
            hdf5_repo.load(sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index)
            is None
        )

        # 他の結果は残っていることを確認
        assert (
            hdf5_repo.load(sample_result2.protein_id, sample_result2.compound_set_id, sample_result2.compound_index)
            is not None
        )

        # HDF5ファイル内部のグループが削除されているか確認 (オプション)
        # グループパスを修正
        group_path1 = (
            f"/results/{sample_result1.protein_id}/{sample_result1.compound_set_id}/{sample_result1.compound_index}"
        )
        group_path2 = (
            f"/results/{sample_result2.protein_id}/{sample_result2.compound_set_id}/{sample_result2.compound_index}"
        )
        with h5py.File(hdf5_repo.hdf5_file_path, "r") as f:
            assert group_path1 not in f
            assert group_path2 in f

    def test_delete_non_existent(self, hdf5_repo: HDF5DockingResultRepository):
        """存在しない結果を削除しようとしてもエラーにならないことを確認する。"""
        try:
            # deleteメソッドの引数を修正
            hdf5_repo.delete("non_existent_protein", "non_existent_set", 0)
        except Exception as e:
            pytest.fail(f"Deleting non-existent entry raised an exception: {e}")

    def test_update(self, hdf5_repo: HDF5DockingResultRepository, sample_result1: DockingResult, tmp_path: Path):
        """updateメソッドがsaveと同じように上書き動作をすることを確認する。"""
        hdf5_repo.save(sample_result1)

        # 更新用のDockingResultを作成
        updated_result_path = tmp_path / "updated_via_update_proteinA_set1_0.sdf"
        updated_result_path.touch()
        updated_result = DockingResult.create(
            result_path=updated_result_path,
            protein_id=sample_result1.protein_id,
            compound_set_id=sample_result1.compound_set_id,
            compound_index=sample_result1.compound_index,
            docking_score=-11.0,
            metadata={"pose_data": "updated via update method"},
        )
        hdf5_repo.update(updated_result)  # updateは内部でsaveを呼ぶ

        # loadメソッドの引数を修正
        loaded_result = hdf5_repo.load(
            sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index
        )
        assert loaded_result is not None
        assert loaded_result.docking_score == -11.0
        assert loaded_result.get_metadata_value("pose_data") == "updated via update method"
        assert loaded_result.result_path == updated_result_path

    def test_get_repository_type(self, hdf5_repo: HDF5DockingResultRepository):
        """get_repository_typeが正しく"hdf5"を返すことを確認する。"""
        assert hdf5_repo.get_repository_type() == "hdf5"

    def test_get_connection_details(self, hdf5_repo: HDF5DockingResultRepository):
        """get_connection_detailsがファイルパスを含む辞書を返すことを確認する。"""
        details = hdf5_repo.get_connection_details()
        assert "file_path" in details
        assert details["file_path"] == str(hdf5_repo.hdf5_file_path)


# --- 並行性テスト用ヘルパー関数 ---
# クラスの外に定義
def _save_task_func(repo_path: Path, result: DockingResult, delay: float = 0):
    """別プロセスで実行される保存タスク。"""
    time.sleep(delay)  # 他のプロセスがロックを取得するのを待つための遅延
    repo = HDF5DockingResultRepository(repo_path)
    repo.save(result)


class TestHDF5DockingResultRepository:
    """HDF5DockingResultRepositoryのテストクラス。"""

    # ... (他のテストメソッドは省略) ...

    # --- 並行性テスト ---
    # 注意: このテストは基本的なロックの動作を確認するものであり、
    #       厳密な並行性問題を完全に検出するものではありません。

    def test_concurrent_save(self, tmp_path: Path, sample_result1: DockingResult, sample_result2: DockingResult):
        """複数のプロセスが同時に保存しようとしてもデータが壊れないことを確認する（基本的なロック確認）。"""
        hdf5_file = tmp_path / "concurrent_test.hdf5"
        repo = HDF5DockingResultRepository(hdf5_file_path=hdf5_file)

        # 異なる結果を別々のプロセスで保存
        # targetをクラス外の関数に変更
        process1 = Process(target=_save_task_func, args=(hdf5_file, sample_result1, 0))
        process2 = Process(target=_save_task_func, args=(hdf5_file, sample_result2, 0.1))  # 少し遅延させる

        process1.start()
        process2.start()

        process1.join(timeout=10)  # タイムアウトを設定
        process2.join(timeout=10)

        assert process1.exitcode == 0, "Process 1 failed"
        assert process2.exitcode == 0, "Process 2 failed"

        # 両方の結果が正しく保存されているか確認
        # loadメソッドの引数を修正
        loaded1 = repo.load(sample_result1.protein_id, sample_result1.compound_set_id, sample_result1.compound_index)
        loaded2 = repo.load(sample_result2.protein_id, sample_result2.compound_set_id, sample_result2.compound_index)

        assert loaded1 is not None
        assert loaded1.docking_score == sample_result1.docking_score
        assert loaded1.get_metadata_value("pose_data") == sample_result1.get_metadata_value("pose_data")
        assert loaded2 is not None
        assert loaded2.docking_score == sample_result2.docking_score
        assert loaded2.get_metadata_value("pose_data") == sample_result2.get_metadata_value("pose_data")

        # ファイルが破損していないか基本的なチェック
        # グループパスを修正
        group_path1 = (
            f"/results/{sample_result1.protein_id}/{sample_result1.compound_set_id}/{sample_result1.compound_index}"
        )
        group_path2 = (
            f"/results/{sample_result2.protein_id}/{sample_result2.compound_set_id}/{sample_result2.compound_index}"
        )
        try:
            with h5py.File(hdf5_file, "r") as f:
                assert group_path1 in f
                assert group_path2 in f
        except Exception as e:
            pytest.fail(f"HDF5 file seems corrupted after concurrent writes: {e}")
