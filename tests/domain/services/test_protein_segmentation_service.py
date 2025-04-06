import os
import shutil
import tempfile
from pathlib import Path
from typing import List
from unittest import mock

import pytest

from docking_automation.domain.services.protein_segmentation_service import (
    ProteinSegmentationService,
)
from docking_automation.molecule.protein import Protein
from docking_automation.molecule.protein_events import ProteinSegmented


class TestProteinSegmentationService:
    """
    ProteinSegmentationServiceのテストクラス。
    """

    @pytest.fixture
    def service(self):
        """
        テスト用のProteinSegmentationServiceインスタンスを提供する。
        """
        return ProteinSegmentationService()

    @pytest.fixture
    def temp_dir(self):
        """
        テスト用の一時ディレクトリを提供する。
        """
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        # テスト後にディレクトリを削除
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def sample_protein(self, temp_dir):
        """
        テスト用のサンプルProteinオブジェクトを提供する。
        """
        # サンプルPDBファイルの作成
        pdb_content = """ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00 16.77           N  
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00 16.57           C  
ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00 16.16           C  
ATOM      4  O   ALA A   1       2.032  -0.994   0.169  1.00 16.78           O  
ATOM      5  CB  ALA A   1      -0.509   0.856   0.727  1.00 16.33           C  
END
"""
        pdb_path = temp_dir / "sample.pdb"
        with open(pdb_path, "w") as f:
            f.write(pdb_content)

        return Protein.create(path=pdb_path)

    @pytest.fixture
    def mock_alpha_cutter_output(self, temp_dir):
        """
        モック用のAlphaCutter出力ファイルを作成する。
        """
        # ドメイン分割されたPDBファイルの作成
        domain1_content = """ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00 16.77           N  
ATOM      2  CA  ALA A   1      -0.001   0.064  -0.491  1.00 16.57           C  
END
"""
        domain2_content = """ATOM      3  C   ALA A   1       1.499  -0.110  -0.491  1.00 16.16           C  
ATOM      4  O   ALA A   1       2.032  -0.994   0.169  1.00 16.78           O  
END
"""

        # サマリーファイルの作成
        summary_content = "protein_id,domain_id,start_residue,end_residue\nsample,domain1,1,2\nsample,domain2,3,4\n"

        # ファイルの作成
        domain1_path = temp_dir / "sample_domain1.pdb"
        domain2_path = temp_dir / "sample_domain2.pdb"
        summary_path = temp_dir / "AFCT-OUT_summary.csv"

        with open(domain1_path, "w") as f:
            f.write(domain1_content)
        with open(domain2_path, "w") as f:
            f.write(domain2_content)
        with open(summary_path, "w") as f:
            f.write(summary_content)

        return [domain1_path, domain2_path, summary_path]

    def test_get_script_path(self, service):
        """
        _get_script_pathメソッドが正しくスクリプトパスを返すことを確認する。
        """
        script_path = service._get_script_path()
        assert script_path.exists()
        assert script_path.name == "alpha_cutter.py"

    def test_build_alpha_cutter_command(self, service, temp_dir):
        """
        _build_alpha_cutter_commandメソッドが正しくコマンドを構築することを確認する。
        """
        script_path = Path("/path/to/alpha_cutter.py")
        input_file_path = temp_dir / "sample.pdb"
        options = {"loop_min": 5, "domain_min": 30, "domain_out": True}

        command = service._build_alpha_cutter_command(script_path, input_file_path, options)

        assert command[0] == "python"
        assert command[1] == str(script_path)
        assert "--loop_min" in command
        assert "5" in command
        assert "--domain_min" in command
        assert "30" in command
        assert "--domain_out" in command

    @mock.patch("docking_automation.domain.services.protein_segmentation_service.subprocess.run")
    def test_run_alpha_cutter(self, mock_run, service, temp_dir):
        """
        _run_alpha_cutterメソッドが正しくsubprocess.runを呼び出すことを確認する。
        """
        command = ["python", "/path/to/alpha_cutter.py", "--domain_out"]

        # モックの設定
        mock_run.return_value.stdout = "AlphaCutter実行完了"
        mock_run.return_value.stderr = ""

        service._run_alpha_cutter(command, temp_dir)

        # subprocess.runが正しく呼び出されたことを確認
        mock_run.assert_called_once_with(command, cwd=str(temp_dir), check=True, capture_output=True, text=True)

    def test_collect_outputs(self, service, temp_dir, mock_alpha_cutter_output):
        """
        _collect_outputsメソッドが正しく出力ファイルを収集することを確認する。
        """
        final_output_dir = temp_dir / "output"
        os.makedirs(final_output_dir, exist_ok=True)

        output_file_paths = service._collect_outputs(temp_dir, final_output_dir)

        # PDBファイルのみが返されることを確認
        assert len(output_file_paths) == 2
        assert all(path.suffix.lower() == ".pdb" for path in output_file_paths)

        # ファイルが移動されたことを確認
        assert all(path.exists() for path in output_file_paths)
        assert all(path.parent == final_output_dir for path in output_file_paths)

    def test_create_protein_objects(self, service, temp_dir, sample_protein):
        """
        _create_protein_objectsメソッドが正しくProteinオブジェクトを作成することを確認する。
        """
        # テスト用のPDBファイル
        file_paths = [temp_dir / "sample_domain1.pdb", temp_dir / "sample_domain2.pdb"]

        # ファイルの作成
        for path in file_paths:
            with open(path, "w") as f:
                f.write("ATOM      1  N   ALA A   1      -0.677  -1.230  -0.491  1.00 16.77           N\nEND\n")

        proteins = service._create_protein_objects(file_paths, sample_protein)

        # 正しい数のProteinオブジェクトが作成されることを確認
        assert len(proteins) == 2

        # IDが正しく設定されていることを確認
        assert proteins[0].id.startswith(sample_protein.id)
        assert proteins[1].id.startswith(sample_protein.id)

        # パスが正しく設定されていることを確認
        assert proteins[0].path == file_paths[0]
        assert proteins[1].path == file_paths[1]

    @mock.patch(
        "docking_automation.domain.services.protein_segmentation_service.ProteinSegmentationService._get_script_path"
    )
    @mock.patch(
        "docking_automation.domain.services.protein_segmentation_service.ProteinSegmentationService._run_alpha_cutter"
    )
    @mock.patch("docking_automation.domain.services.protein_segmentation_service.shutil.copy")
    @mock.patch("docking_automation.domain.services.protein_segmentation_service.shutil.move")
    def test_segment(
        self,
        mock_move,
        mock_copy,
        mock_run,
        mock_get_script_path,
        service,
        temp_dir,
        sample_protein,
        mock_alpha_cutter_output,
    ):
        """
        segmentメソッドが正しく動作することを確認する。
        """
        # モックの設定
        mock_get_script_path.return_value = Path("/path/to/alpha_cutter.py")

        # 出力ディレクトリ
        final_output_dir = temp_dir / "output"
        os.makedirs(final_output_dir, exist_ok=True)

        # テスト対象メソッドの呼び出し
        with mock.patch("docking_automation.domain.services.protein_segmentation_service.glob.glob") as mock_glob:
            # globの戻り値を設定
            mock_glob.side_effect = lambda pattern: [
                str(path) for path in mock_alpha_cutter_output if path.match(pattern.replace("*", ""))
            ]

            # tempfile.TemporaryDirectoryのモック
            with mock.patch(
                "docking_automation.domain.services.protein_segmentation_service.tempfile.TemporaryDirectory"
            ) as mock_temp_dir:
                mock_temp_dir.return_value.__enter__.return_value = str(temp_dir)

                # Path.existsのモックを先に適用
                with mock.patch("pathlib.Path.exists") as mock_exists:
                    mock_exists.return_value = True

                    # Protein.createのモック
                    with mock.patch("docking_automation.molecule.protein.Protein.create") as mock_create:
                        # モックProteinオブジェクトを返す
                        # 複数のモックProteinオブジェクトを作成
                        mock_proteins = []
                        for i in range(2):  # 2つのドメインを想定
                            mock_protein = mock.MagicMock(spec=Protein)
                            mock_protein.id = f"mocked_protein_id_{i}"
                            mock_protein.path = Path(f"/mocked/path_{i}.pdb")
                            mock_proteins.append(mock_protein)

                        # 最初のProteinオブジェクトにProteinSegmentedイベントを設定
                        event = ProteinSegmented(
                            original_protein_id=sample_protein.id, segmented_protein_ids=[p.id for p in mock_proteins]
                        )
                        mock_proteins[0].release_domain_events.return_value = [event]

                        # 他のProteinオブジェクトは空のイベントリストを返す
                        for i in range(1, len(mock_proteins)):
                            mock_proteins[i].release_domain_events.return_value = []

                        # create呼び出しごとに異なるモックを返す
                        mock_create.side_effect = mock_proteins
                        mock_create.return_value = mock_protein

                        options = {"domain_out": True, "loop_min": 5}
                        result = service.segment(sample_protein, options, final_output_dir)

        # 結果の検証
        assert len(result) > 0

        # ドメインイベントが発行されていることを確認
        events = result[0].release_domain_events()
        assert any(isinstance(event, ProteinSegmented) for event in events)

        # ProteinSegmentedイベントの内容を確認
        for event in events:
            if isinstance(event, ProteinSegmented):
                assert event.original_protein_id == sample_protein.id
                assert len(event.segmented_protein_ids) > 0
