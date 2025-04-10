"""
HDF5リポジトリの上書きモードと追記モードの使用例。

このサンプルでは、HDF5DockingResultRepositoryの上書きモードと追記モードの違いを示します。
"""

import os
from pathlib import Path

from docking_automation.docking.docking_result import DockingResult
from docking_automation.infrastructure.repositories.docking_result_repository_factory import (
    DockingResultRepositoryFactory,
    RepositoryType,
)
from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import (
    HDF5DockingResultRepository,
)


def main():
    """上書きモードと追記モードの違いを示すメイン関数。"""
    # 一時ディレクトリを作成
    output_dir = Path("output/hdf5_example")
    output_dir.mkdir(parents=True, exist_ok=True)

    # サンプルのドッキング結果を作成
    result_file1 = output_dir / "protein1_compound1.sdf"
    result_file1.touch()
    result1 = DockingResult.create(
        result_path=result_file1,
        protein_id="protein1",
        compound_set_id="compound_set1",
        compound_index=0,
        docking_score=-9.5,
        metadata={"pose_data": "Original pose data"},
    )

    # 更新用のドッキング結果を作成
    result_file2 = output_dir / "protein1_compound1_updated.sdf"
    result_file2.touch()
    result2 = DockingResult.create(
        result_path=result_file2,
        protein_id="protein1",
        compound_set_id="compound_set1",
        compound_index=0,  # 同じキー
        docking_score=-10.5,  # 異なるスコア
        metadata={"pose_data": "Updated pose data"},  # 異なるメタデータ
    )

    # 1. 上書きモードの例
    print("===== 上書きモードの例 =====")
    overwrite_file = output_dir / "overwrite_example.hdf5"
    if overwrite_file.exists():
        os.remove(overwrite_file)

    # ファクトリを使用してリポジトリを作成（上書きモード）
    overwrite_repo = DockingResultRepositoryFactory.create(
        RepositoryType.HDF5,
        base_directory=output_dir,
        config={"mode": "overwrite"},
    )
    # または直接インスタンス化
    # overwrite_repo = HDF5DockingResultRepository(hdf5_file_path=overwrite_file, mode="overwrite")

    # 最初の結果を保存
    print("最初の結果を保存します...")
    overwrite_repo.save(result1)

    # 保存された結果を読み込み
    loaded_result = overwrite_repo.load(result1.protein_id, result1.compound_set_id, result1.compound_index)
    print(f"読み込んだ結果: スコア={loaded_result.docking_score}, メタデータ={loaded_result.metadata}")

    # 更新された結果を保存（上書きモード）
    print("更新された結果を保存します（上書きモード）...")
    overwrite_repo.save(result2)

    # 更新された結果を読み込み
    loaded_result = overwrite_repo.load(result1.protein_id, result1.compound_set_id, result1.compound_index)
    print(f"読み込んだ結果: スコア={loaded_result.docking_score}, メタデータ={loaded_result.metadata}")
    print("上書きモードでは、既存の結果が上書きされます。")

    # 2. 追記モードの例
    print("\n===== 追記モードの例 =====")
    append_file = output_dir / "append_example.hdf5"
    if append_file.exists():
        os.remove(append_file)

    # ファクトリを使用してリポジトリを作成（追記モード）
    append_repo = DockingResultRepositoryFactory.create(
        RepositoryType.HDF5,
        base_directory=output_dir,
        config={"mode": "append"},
    )
    # または直接インスタンス化
    # append_repo = HDF5DockingResultRepository(hdf5_file_path=append_file, mode="append")

    # 最初の結果を保存
    print("最初の結果を保存します...")
    append_repo.save(result1)

    # 保存された結果を読み込み
    loaded_result = append_repo.load(result1.protein_id, result1.compound_set_id, result1.compound_index)
    print(f"読み込んだ結果: スコア={loaded_result.docking_score}, メタデータ={loaded_result.metadata}")

    # 更新された結果を保存（追記モード）
    print("更新された結果を保存します（追記モード）...")
    append_repo.save(result2)

    # 結果を読み込み（元のデータが保持されている）
    loaded_result = append_repo.load(result1.protein_id, result1.compound_set_id, result1.compound_index)
    print(f"読み込んだ結果: スコア={loaded_result.docking_score}, メタデータ={loaded_result.metadata}")
    print("追記モードでは、既存の結果が保持されます。")

    # 3. updateメソッドの例（モードに関わらず常に上書き）
    print("\n===== updateメソッドの例 =====")
    # 追記モードのリポジトリでupdateメソッドを使用
    print("追記モードのリポジトリでupdateメソッドを使用...")
    append_repo.update(result2)  # updateは常に上書きモードで動作

    # 結果を読み込み（updateメソッドにより上書きされている）
    loaded_result = append_repo.load(result1.protein_id, result1.compound_set_id, result1.compound_index)
    print(f"読み込んだ結果: スコア={loaded_result.docking_score}, メタデータ={loaded_result.metadata}")
    print("updateメソッドは、モード設定に関わらず常に上書きモードで動作します。")


if __name__ == "__main__":
    main()
