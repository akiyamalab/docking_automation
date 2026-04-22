from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from docking_automation.docking.fpocket_grid_box_predictor import FpocketGridBoxPredictor


def _make_protein_mock(path: Path) -> MagicMock:
    """Proteinモックを作成する"""
    protein = MagicMock()
    protein.path = path
    protein.id = path.stem
    return protein


def _build_fake_fpocket_output(protein_path: Path) -> None:
    """fpocketが生成するディレクトリ・ファイルを模倣して作成する"""
    output_name = protein_path.stem
    fpocket_output_dir = protein_path.parent / f"{output_name}_out"
    fpocket_output_dir.mkdir(parents=True, exist_ok=True)
    output_pdb = fpocket_output_dir / f"{output_name}_out.pdb"
    # STP残基を1つ含む最小PDBファイル
    output_pdb.write_text(
        "HETATM    1  O   STP     1      10.000  20.000  30.000  0.00  0.00          O\n"
        "END\n"
    )


class TestFpocketGridBoxPredictor:
    def test_predict_with_absolute_path(self, tmp_path):
        """絶対パスを持つProteinでfpocketが呼び出せること（相対パスエラーが出ないこと）"""
        # 絶対パスのダミーPDBファイルを作成
        protein_pdb = tmp_path / "protein.pdb"
        protein_pdb.write_text("ATOM\n")

        protein = _make_protein_mock(protein_pdb)
        predictor = FpocketGridBoxPredictor()

        # subprocess.run をモック: 成功扱いにし、fpocket出力を模倣
        def fake_run(cmd, **kwargs):
            _build_fake_fpocket_output(protein_pdb)
            return MagicMock(returncode=0, stdout="", stderr="")

        with patch("docking_automation.docking.fpocket_grid_box_predictor.subprocess.run", side_effect=fake_run):
            try:
                result = predictor.predict(protein, pocket_rank=1)
                # 正常系: GridBoxが返るはず
                assert result is not None
            except ValueError as e:
                # fpocket実行エラーは許容するが、relative_to由来のエラーは禁止
                assert "is not relative to" not in str(e), (
                    f"relative_to() バグが再発しました: {e}"
                )

    def test_predict_with_deeply_nested_absolute_path(self, tmp_path):
        """深いネストの絶対パスでもrelative_to ValueErrorが出ないこと"""
        nested = tmp_path / "a" / "b" / "c"
        nested.mkdir(parents=True)
        protein_pdb = nested / "deep_protein.pdb"
        protein_pdb.write_text("ATOM\n")

        protein = _make_protein_mock(protein_pdb)
        predictor = FpocketGridBoxPredictor()

        def fake_run(cmd, **kwargs):
            _build_fake_fpocket_output(protein_pdb)
            return MagicMock(returncode=0, stdout="", stderr="")

        with patch("docking_automation.docking.fpocket_grid_box_predictor.subprocess.run", side_effect=fake_run):
            try:
                result = predictor.predict(protein, pocket_rank=1)
                assert result is not None
            except ValueError as e:
                assert "is not relative to" not in str(e), (
                    f"relative_to() バグが再発しました: {e}"
                )
