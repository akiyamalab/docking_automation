from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from docking_automation.docking.docking_parameters import (
    CommonDockingParameters,
    DockingParameters,
    UniDockParameters,
)
from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from docking_automation.docking.unidock_docking import UniDockDocking

VINA_RESULT_LINE = "REMARK VINA RESULT:   -8.5      0.000      0.000\n"


@pytest.fixture
def grid_box():
    mock = MagicMock()
    mock.center = [10.0, 20.0, 30.0]
    mock.size = [25.0, 25.0, 25.0]
    return mock


@pytest.fixture
def preprocessed_protein(tmp_path):
    receptor = tmp_path / "protein.pdbqt"
    receptor.write_text("fake receptor content")
    protein = MagicMock(spec=PreprocessedProtein)
    protein.file_path = receptor
    protein.content_hash = "protein_hash_abc"
    return protein


@pytest.fixture
def ligand_paths(tmp_path):
    paths = []
    for i in range(2):
        p = tmp_path / f"compound_{i}.pdbqt"
        p.write_text(f"fake ligand {i}")
        paths.append(p)
    return paths


@pytest.fixture
def preprocessed_compound_set(ligand_paths):
    cs = MagicMock(spec=PreprocessedCompoundSet)
    cs.file_paths = ligand_paths
    cs.get_compound_hash.side_effect = lambda idx: f"hash_{idx}"
    return cs


def _make_fake_run(ligand_count=None, score_line=VINA_RESULT_LINE):
    """fake subprocess.run factory: out_dir に {stem}_out.pdbqt を生成する"""
    def fake_run(cmd, **kwargs):
        out_dir = Path(cmd[cmd.index("--dir") + 1])
        ligand_idx_path = Path(cmd[cmd.index("--ligand_index") + 1])
        paths = ligand_idx_path.read_text().splitlines()
        count = ligand_count if ligand_count is not None else len(paths)
        for ligand_path in paths[:count]:
            stem = Path(ligand_path).stem
            (out_dir / f"{stem}_out.pdbqt").write_text(score_line)
        return MagicMock(returncode=0, stdout="", stderr="")
    return fake_run


class TestUniDockParameters:
    def test_unidock_parameters_defaults(self):
        params = UniDockParameters()
        assert params.seed == 1
        assert params.search_mode == "balance"
        assert params.scoring == "vina"
        assert params.num_modes == 1
        assert params.rescue_mode is False
        assert params.rescue_search_mode == "detail"


class TestUniDockDocking:
    def test_build_cli_command(self, preprocessed_protein, preprocessed_compound_set, grid_box, tmp_path):
        docking = UniDockDocking()
        ligand_index = tmp_path / "ligands.txt"
        ligand_index.write_text("\n".join(str(p) for p in preprocessed_compound_set.file_paths))
        out_dir = tmp_path / "output"
        out_dir.mkdir()
        params = UniDockParameters()

        cmd = docking._build_cli_command(
            preprocessed_protein.file_path,
            ligand_index,
            grid_box,
            out_dir,
            params,
        )

        # 必要なフラグが全て含まれていることを確認
        assert "--receptor" in cmd and str(preprocessed_protein.file_path) in cmd
        assert "--ligand_index" in cmd and str(ligand_index) in cmd
        assert "--center_x" in cmd and str(grid_box.center[0]) in cmd
        assert "--center_y" in cmd and str(grid_box.center[1]) in cmd
        assert "--center_z" in cmd and str(grid_box.center[2]) in cmd
        assert "--size_x" in cmd and str(grid_box.size[0]) in cmd
        assert "--size_y" in cmd and str(grid_box.size[1]) in cmd
        assert "--size_z" in cmd and str(grid_box.size[2]) in cmd
        assert "--scoring" in cmd and params.scoring in cmd
        assert "--search_mode" in cmd and params.search_mode in cmd
        assert "--num_modes" in cmd and str(params.num_modes) in cmd
        assert "--seed" in cmd and str(params.seed) in cmd

    def test_run_docking_success(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        docking = UniDockDocking()

        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=_make_fake_run()):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None

                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(common=common, specific=UniDockParameters()))

        assert len(results) == 2
        for result in results:
            assert isinstance(result, DockingResult)
            assert result.docking_score == pytest.approx(-8.5)

    def test_run_docking_failed_ligand(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """一部リガンドの出力ファイルが存在しない場合、失敗リガンドはresultsから除外される"""
        docking = UniDockDocking()

        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=_make_fake_run(ligand_count=1)):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None

                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(common=common, specific=UniDockParameters()))

        assert len(results) == 1
        assert results[0].docking_score == pytest.approx(-8.5)

    def test_penalty_score_filtered(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """ペナルティスコア(999999)はフィルタされ結果に含まれない"""
        penalty_line = "REMARK VINA RESULT:   999999      0.000      0.000\n"
        docking = UniDockDocking()

        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=_make_fake_run(score_line=penalty_line)):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(common=common, specific=UniDockParameters()))

        assert len(results) == 0

    def test_valid_score_passes(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """正常スコア(-7.5)はフィルタされず結果に含まれる"""
        valid_line = "REMARK VINA RESULT:   -7.5      0.000      0.000\n"
        docking = UniDockDocking()

        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=_make_fake_run(score_line=valid_line)):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(common=common, specific=UniDockParameters()))

        assert len(results) == 2
        for result in results:
            assert result.docking_score == pytest.approx(-7.5)

    def test_threshold_boundary(self, preprocessed_protein, grid_box, tmp_path):
        """境界値: score=5.0 → フィルタ(None)、score=4.9 → 通過"""
        ligand_paths = []
        score_lines = ["REMARK VINA RESULT:   5.0      0.000      0.000\n",
                       "REMARK VINA RESULT:   4.9      0.000      0.000\n"]
        for i, line in enumerate(score_lines):
            p = tmp_path / f"compound_{i}.pdbqt"
            p.write_text(f"fake ligand {i}")
            ligand_paths.append(p)

        cs = MagicMock(spec=PreprocessedCompoundSet)
        cs.file_paths = ligand_paths
        cs.get_compound_hash.side_effect = lambda idx: f"hash_{idx}"

        def fake_run_multi(cmd, **kwargs):
            out_dir = Path(cmd[cmd.index("--dir") + 1])
            ligand_idx_path = Path(cmd[cmd.index("--ligand_index") + 1])
            paths = ligand_idx_path.read_text().splitlines()
            for j, ligand_path in enumerate(paths):
                stem = Path(ligand_path).stem
                (out_dir / f"{stem}_out.pdbqt").write_text(score_lines[j])
            return MagicMock(returncode=0, stdout="", stderr="")

        docking = UniDockDocking()
        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=fake_run_multi):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=cs,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(common=common, specific=UniDockParameters()))

        assert len(results) == 1
        assert results[0].docking_score == pytest.approx(4.9)

    def test_rescue_mode_retries_failed_ligands(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """rescue_mode=True のとき、failed リガンドが rescue_search_mode で再実行される"""
        call_count = [0]

        def fake_run_with_rescue(cmd, **kwargs):
            call_count[0] += 1
            out_dir = Path(cmd[cmd.index("--dir") + 1])
            ligand_idx_path = Path(cmd[cmd.index("--ligand_index") + 1])
            paths = ligand_idx_path.read_text().splitlines()
            if call_count[0] == 1:
                # 1回目: compound_0 のみ成功 (compound_1 は出力ファイルなし)
                stem = Path(paths[0]).stem
                (out_dir / f"{stem}_out.pdbqt").write_text(VINA_RESULT_LINE)
            else:
                # rescue: 渡された全リガンドが成功
                for p in paths:
                    stem = Path(p).stem
                    (out_dir / f"{stem}_out.pdbqt").write_text(VINA_RESULT_LINE)
            return MagicMock(returncode=0, stdout="", stderr="")

        docking = UniDockDocking()
        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=fake_run_with_rescue):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(
                    common=common,
                    specific=UniDockParameters(rescue_mode=True, rescue_search_mode="detail"),
                ))

        assert len(results) == 2
        assert call_count[0] == 2

    def test_rescue_mode_false_no_retry(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """rescue_mode=False のとき、failed リガンドは再実行されない"""
        call_count = [0]

        def fake_run_partial(cmd, **kwargs):
            call_count[0] += 1
            out_dir = Path(cmd[cmd.index("--dir") + 1])
            ligand_idx_path = Path(cmd[cmd.index("--ligand_index") + 1])
            paths = ligand_idx_path.read_text().splitlines()
            # compound_0 のみ成功
            stem = Path(paths[0]).stem
            (out_dir / f"{stem}_out.pdbqt").write_text(VINA_RESULT_LINE)
            return MagicMock(returncode=0, stdout="", stderr="")

        docking = UniDockDocking()
        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=fake_run_partial):
            with patch.object(docking, "converter") as mock_converter:
                mock_converter.pdbqt_to_sdf.return_value = None
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                results = docking.dock(DockingParameters(
                    common=common,
                    specific=UniDockParameters(rescue_mode=False),
                ))

        assert len(results) == 1
        assert call_count[0] == 1

    def test_run_docking_uses_seed(self, preprocessed_protein, preprocessed_compound_set, grid_box):
        """run_docking() 呼び出し時に seed=1 が CLI に渡されることを確認"""
        docking = UniDockDocking()
        captured_cmd: list = []

        def spy_run(cmd, **kwargs):
            captured_cmd.extend(cmd)
            return MagicMock(returncode=0, stdout="", stderr="")

        with patch("docking_automation.docking.unidock_docking.subprocess.run", side_effect=spy_run):
            with patch.object(docking, "converter"):
                common = CommonDockingParameters(
                    protein=preprocessed_protein,
                    compound_set=preprocessed_compound_set,
                    grid_box=grid_box,
                )
                docking.dock(DockingParameters(common=common, specific=UniDockParameters(seed=1)))

        assert "--seed" in captured_cmd
        assert captured_cmd[captured_cmd.index("--seed") + 1] == "1"
