"""
ドッキング結果エクスポートサービスのテスト。
"""

import tempfile
from pathlib import Path

import pytest

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.domain.services.docking_result_export_service import (
    export_docking_results_for_protein,
)


class TestDockingResultExportService:
    """ドッキング結果エクスポートサービスのテストクラス。"""

    def test_export_docking_results_for_protein(self, mocker):
        """タンパク質に対するドッキング結果をエクスポートするテスト。"""
        # モックリポジトリを作成
        mock_repository = mocker.Mock()

        # テスト用のドッキング結果を作成
        result1 = DockingResult(
            result_path=Path("test1.pdbqt"),
            protein_id="protein1",
            compound_set_id="set1",
            compound_index=0,
            docking_score=-8.5,
            protein_content_hash="protein_hash_1",
            compoundset_content_hash="compound_hash_1",
            metadata={"ligand_efficiency": 0.5},
        )

        result2 = DockingResult(
            result_path=Path("test2.pdbqt"),
            protein_id="protein1",
            compound_set_id="set1",
            compound_index=1,
            docking_score=-7.2,
            protein_content_hash="protein_hash_1",
            compoundset_content_hash="compound_hash_1",
            metadata={"ligand_efficiency": 0.4},
        )

        # DockingResultCollectionを作成
        collection = DockingResultCollection([result1, result2])

        # find_by_proteinメソッドのモックを設定
        mock_repository.find_by_protein.return_value = collection

        # export_to_sdfメソッドのモックを設定
        mocker.patch.object(collection, "export_to_sdf")

        # 一時ファイルのパスを作成
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
            temp_path = Path(temp_file.name)

        try:
            # ドッキング結果をエクスポート
            result_path = export_docking_results_for_protein(
                repository=mock_repository, protein_id="protein1", output_path=temp_path
            )

            # find_by_proteinメソッドが呼ばれたことを確認
            mock_repository.find_by_protein.assert_called_once_with("protein1")

            # export_to_sdfメソッドが呼ばれたことを確認
            collection.export_to_sdf.assert_called_once_with(temp_path)

            # 戻り値が正しいことを確認
            assert result_path == temp_path
        finally:
            # 一時ファイルを削除
            if temp_path.exists():
                temp_path.unlink()

    def test_export_docking_results_for_protein_default_path(self, mocker):
        """デフォルトのパスを使用してエクスポートするテスト。"""
        # モックリポジトリを作成
        mock_repository = mocker.Mock()

        # テスト用のドッキング結果を作成
        result = DockingResult(
            result_path=Path("test.pdbqt"),
            protein_id="protein1",
            compound_set_id="set1",
            compound_index=0,
            docking_score=-8.5,
            protein_content_hash="protein_hash_1",
            compoundset_content_hash="compound_hash_1",
        )

        # DockingResultCollectionを作成
        collection = DockingResultCollection([result])

        # find_by_proteinメソッドのモックを設定
        mock_repository.find_by_protein.return_value = collection

        # export_to_sdfメソッドのモックを設定
        mocker.patch.object(collection, "export_to_sdf")

        # ドッキング結果をエクスポート（出力パスを指定しない）
        result_path = export_docking_results_for_protein(repository=mock_repository, protein_id="protein1")

        # find_by_proteinメソッドが呼ばれたことを確認
        mock_repository.find_by_protein.assert_called_once_with("protein1")

        # export_to_sdfメソッドが呼ばれたことを確認
        expected_path = Path("protein1_docking_results.sdf")
        collection.export_to_sdf.assert_called_once()
        actual_path = collection.export_to_sdf.call_args[0][0]
        assert actual_path.name == expected_path.name

        # 戻り値が正しいことを確認
        assert result_path.name == expected_path.name

        # 一時ファイルを削除（もし作成されていれば）
        if result_path.exists():
            result_path.unlink()

    def test_export_docking_results_for_protein_empty_results(self, mocker):
        """結果が空の場合のテスト。"""
        # モックリポジトリを作成
        mock_repository = mocker.Mock()

        # 空のDockingResultCollectionを作成
        collection = DockingResultCollection()

        # find_by_proteinメソッドのモックを設定
        mock_repository.find_by_protein.return_value = collection

        # 結果が空の場合はValueErrorが発生することを確認
        with pytest.raises(ValueError) as excinfo:
            export_docking_results_for_protein(repository=mock_repository, protein_id="protein1")

        # エラーメッセージを確認
        assert "タンパク質ID 'protein1' に対するドッキング結果が見つかりません" in str(excinfo.value)
