from pathlib import Path

import pytest

from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
    DockingResultRepositoryFactory,
    RepositoryType,
)
from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
    HDF5DockingResultRepository,
)


class TestDockingResultRepositoryFactory:
    """DockingResultRepositoryFactoryのテストクラス。"""

    def test_create_hdf5_repository_with_mode(self, tmp_path: Path):
        """HDF5リポジトリの作成時にモードを指定できることを確認する。"""
        # 上書きモード（デフォルト）
        repo1 = DockingResultRepositoryFactory.create(
            RepositoryType.HDF5,
            base_directory=tmp_path / "repo1",
        )
        assert isinstance(repo1, HDF5DockingResultRepository)
        assert repo1.mode == "overwrite"  # デフォルトは上書きモード

        # 上書きモード（明示的に指定）
        repo2 = DockingResultRepositoryFactory.create(
            RepositoryType.HDF5,
            base_directory=tmp_path / "repo2",
            config={"mode": "overwrite"},
        )
        assert isinstance(repo2, HDF5DockingResultRepository)
        assert repo2.mode == "overwrite"

        # 追記モード
        repo3 = DockingResultRepositoryFactory.create(
            RepositoryType.HDF5,
            base_directory=tmp_path / "repo3",
            config={"mode": "append"},
        )
        assert isinstance(repo3, HDF5DockingResultRepository)
        assert repo3.mode == "append"

    def test_create_default_with_mode(self):
        """create_defaultメソッドでモードを指定できることを確認する。"""
        # デフォルトモード（上書き）
        repo1 = DockingResultRepositoryFactory.create_default()
        assert isinstance(repo1, HDF5DockingResultRepository)
        assert repo1.mode == "overwrite"

        # 上書きモード（明示的に指定）
        repo2 = DockingResultRepositoryFactory.create_default(mode="overwrite")
        assert isinstance(repo2, HDF5DockingResultRepository)
        assert repo2.mode == "overwrite"

        # 追記モード
        repo3 = DockingResultRepositoryFactory.create_default(mode="append")
        assert isinstance(repo3, HDF5DockingResultRepository)
        assert repo3.mode == "append"
