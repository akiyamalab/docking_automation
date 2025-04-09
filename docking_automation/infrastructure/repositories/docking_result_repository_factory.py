"""
ドッキング結果リポジトリファクトリモジュール。

このモジュールは、設定に基づいて適切なリポジトリ実装を作成するためのファクトリクラスを提供します。
"""

from enum import Enum, auto
from pathlib import Path
from typing import Dict, Optional, Union

from docking_automation.infrastructure.repositories.docking_result_repository import (
    DockingResultRepositoryABC,
)
from docking_automation.infrastructure.repositories.file_docking_result_repository import (
    FileDockingResultRepository,
)
from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
    HDF5DockingResultRepository,
)


class RepositoryType(Enum):
    """
    リポジトリの種類を表す列挙型。
    """

    FILE = auto()
    HDF5 = auto()


class DockingResultRepositoryFactory:
    """
    ドッキング結果リポジトリを作成するファクトリクラス。

    設定に基づいて適切なリポジトリ実装を作成します。
    """

    @staticmethod
    def create(
        repository_type: RepositoryType,
        base_directory: Optional[Path] = None,
        config: Optional[Dict[str, Union[str, int, float, bool]]] = None,
    ) -> DockingResultRepositoryABC:
        """
        ドッキング結果リポジトリを作成する。

        Args:
            repository_type: リポジトリの種類
            base_directory: リポジトリのベースディレクトリ（指定しない場合はデフォルトのパスを使用）
            config: リポジトリの設定（指定しない場合はデフォルトの設定を使用）

        Returns:
            作成されたドッキング結果リポジトリ

        Raises:
            ValueError: 未対応のリポジトリ種類が指定された場合
        """
        if base_directory is None:
            # デフォルトのパスを使用
            base_directory = Path.home() / ".docking_automation" / "repository"

        config = config or {}

        if repository_type == RepositoryType.FILE:
            return FileDockingResultRepository.create(base_directory=base_directory)
        elif repository_type == RepositoryType.HDF5:
            hdf5_file_path = base_directory / "docking_results.hdf5"
            # ディレクトリが存在しない場合は作成
            hdf5_file_path.parent.mkdir(parents=True, exist_ok=True)
            return HDF5DockingResultRepository(hdf5_file_path=hdf5_file_path)
        else:
            raise ValueError(f"未対応のリポジトリ種類です: {repository_type}")

    @staticmethod
    def create_default() -> DockingResultRepositoryABC:
        """
        デフォルト設定でドッキング結果リポジトリを作成する。

        Returns:
            作成されたドッキング結果リポジトリ
        """
        # デフォルトをHDF5に変更
        return DockingResultRepositoryFactory.create(RepositoryType.HDF5)
