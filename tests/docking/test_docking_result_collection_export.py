"""
DockingResultCollectionのエクスポート機能のテスト。
"""

import tempfile
from pathlib import Path

import pytest
from rdkit import Chem

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection


class TestDockingResultCollectionExport:
    """DockingResultCollectionのエクスポート機能のテストクラス。"""

    def test_export_to_sdf_empty_collection(self):
        """空のコレクションをエクスポートするテスト。"""
        collection = DockingResultCollection()

        # 一時ファイルに出力
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
            temp_path = Path(temp_file.name)

        try:
            # SDFファイルにエクスポート
            collection.export_to_sdf(temp_path)

            # 出力されたSDFファイルを検証
            suppl = Chem.SDMolSupplier(str(temp_path))
            mols = [mol for mol in suppl if mol is not None]

            # 空のコレクションなので、分子は0個
            assert len(mols) == 0
        finally:
            # 一時ファイルを削除
            if temp_path.exists():
                temp_path.unlink()

    def test_export_to_sdf_with_mock_results(self, mocker):
        """モックを使用してエクスポート機能をテスト。"""
        # MoleculeConverterのモックを作成
        mock_converter = mocker.patch("docking_automation.converters.molecule_converter.MoleculeConverter")
        mock_instance = mock_converter.return_value

        # pdbqt_to_rdkitメソッドのモックを設定
        mock_mol = Chem.MolFromSmiles("CC")  # 簡単な分子を作成
        mock_instance.pdbqt_to_rdkit.return_value = [mock_mol]

        # テスト用のドッキング結果を作成
        result1 = DockingResult(
            result_path=Path("test1.pdbqt"),
            protein_id="protein1",
            compound_set_id="set1",
            compound_index=0,
            docking_score=-8.5,
            metadata={"ligand_efficiency": 0.5},
        )

        result2 = DockingResult(
            result_path=Path("test2.pdbqt"),
            protein_id="protein1",
            compound_set_id="set1",
            compound_index=1,
            docking_score=-7.2,
            metadata={"ligand_efficiency": 0.4},
        )

        # DockingResultCollectionを作成
        collection = DockingResultCollection([result1, result2])

        # 一時ファイルに出力
        with tempfile.NamedTemporaryFile(suffix=".sdf", delete=False) as temp_file:
            temp_path = Path(temp_file.name)

        try:
            # SDFファイルにエクスポート
            collection.export_to_sdf(temp_path)

            # pdbqt_to_rdkitメソッドが2回呼ばれたことを確認
            assert mock_instance.pdbqt_to_rdkit.call_count == 2

            # 出力されたSDFファイルを検証
            suppl = Chem.SDMolSupplier(str(temp_path))
            mols = [mol for mol in suppl if mol is not None]

            # 2つの結果から各1つの分子が出力されるので、合計2つの分子
            assert len(mols) == 2

            # 各分子のプロパティを検証
            for mol in mols:
                assert mol.HasProp("protein_id")
                assert mol.HasProp("compound_set_id")
                assert mol.HasProp("compound_index")
                assert mol.HasProp("docking_score")
                assert mol.HasProp("ligand_efficiency")

                # プロパティの値を検証
                assert mol.GetProp("protein_id") == "protein1"
                assert mol.GetProp("compound_set_id") == "set1"
                assert mol.GetProp("compound_index") in ["0", "1"]
                assert float(mol.GetProp("docking_score")) in [-8.5, -7.2]
        finally:
            # 一時ファイルを削除
            if temp_path.exists():
                temp_path.unlink()
