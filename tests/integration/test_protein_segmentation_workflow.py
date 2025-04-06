import os
import shutil
import tempfile
from pathlib import Path
from typing import List

import pytest

from docking_automation.domain.services.protein_segmentation_service import (
    ProteinSegmentationService,
)
from docking_automation.molecule.protein import Protein
from docking_automation.molecule.protein_events import ProteinSegmented


class TestProteinSegmentationWorkflow:
    """
    AlphaCutterを使用したタンパク質セグメンテーションワークフローの統合テスト。

    実際のAlphaFold由来のPDBファイルを使用して、ProteinSegmentationServiceを通じて
    AlphaCutterを実行し、タンパク質のセグメンテーション（ドメイン分割・不要残基削除）が
    正しく行われることを確認します。
    """

    @pytest.fixture
    def service(self):
        """
        テスト用のProteinSegmentationServiceインスタンスを提供する。
        """
        return ProteinSegmentationService()

    @pytest.fixture
    def temp_output_dir(self):
        """
        テスト用の一時出力ディレクトリを提供する。
        """
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        # テスト後にディレクトリを削除
        shutil.rmtree(temp_dir)

    @pytest.fixture
    def test_proteins(self):
        """
        テスト用のタンパク質ファイルのパスリストを提供する。
        """
        base_dir = Path(__file__).parent / "test_data" / "alphacutter"
        protein_files = [
            base_dir / "AF-A0A0A2JW93-F1-model_v4.pdb",
            base_dir / "AF-A0A0G2KIZ8-F1-model_v4.pdb",
            base_dir / "AF-A0A0H2UNG0-F1-model_v4.pdb",
        ]

        # ファイルの存在確認
        for file_path in protein_files:
            assert file_path.exists(), f"テストファイルが見つかりません: {file_path}"

        return protein_files

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_index", [0, 1, 2])
    def test_protein_segmentation(self, service, temp_output_dir, test_proteins, protein_index):
        """
        タンパク質セグメンテーションの統合テスト。

        AlphaCutterを使用して、タンパク質のセグメンテーション（ドメイン分割・不要残基削除）が
        正しく行われることを確認します。

        Args:
            service: ProteinSegmentationServiceインスタンス
            temp_output_dir: 一時出力ディレクトリ
            test_proteins: テスト用のタンパク質ファイルのパスリスト
            protein_index: テスト対象のタンパク質のインデックス
        """
        # テスト対象のタンパク質ファイル
        protein_path = test_proteins[protein_index]
        protein_id = protein_path.stem

        # Proteinオブジェクトの作成
        protein = Protein.create(path=protein_path, id=protein_id)

        # AlphaCutterのオプション設定（GitHubの実行例に基づく）
        options = {
            "loop_min": 20,
            "helix_min": 30,
            "fragment_min": 5,
            "domain_min": 50,
            "pLDDT_min": 0,
            "local_contact_range": 5,
            "domain_out": True,
            "single_out": True,
        }

        # タンパク質セグメンテーションの実行
        segmented_proteins = service.segment(protein, options, temp_output_dir)

        # 結果の検証
        # 1. 少なくとも1つのセグメントが生成されていること
        assert len(segmented_proteins) > 0, "セグメント化されたタンパク質が生成されていません"

        # 2. 各セグメントのファイルが存在すること
        for segmented_protein in segmented_proteins:
            assert (
                segmented_protein.path.exists()
            ), f"セグメント化されたタンパク質ファイルが見つかりません: {segmented_protein.path}"

        # 3. ドメインイベントが発行されていること
        events = segmented_proteins[0].release_domain_events()
        protein_segmented_events = [e for e in events if isinstance(e, ProteinSegmented)]

        assert len(protein_segmented_events) > 0, "ProteinSegmentedイベントが発行されていません"

        # 4. イベントの内容が正しいこと
        event = protein_segmented_events[0]
        assert event.original_protein_id == protein.id, "イベントの元タンパク質IDが一致しません"
        assert len(event.segmented_protein_ids) == len(segmented_proteins), "イベントのセグメント数が一致しません"

        # 5. サマリーファイルが生成されていること
        summary_file = temp_output_dir / "AFCT-OUT_summary.csv"
        assert summary_file.exists(), "サマリーファイルが生成されていません"

    @pytest.mark.integration
    def test_protein_segmentation_with_invalid_options(self, service, temp_output_dir, test_proteins):
        """
        無効なオプションを指定した場合のテスト。

        Args:
            service: ProteinSegmentationServiceインスタンス
            temp_output_dir: 一時出力ディレクトリ
            test_proteins: テスト用のタンパク質ファイルのパスリスト
        """
        # テスト対象のタンパク質ファイル
        protein_path = test_proteins[0]
        protein_id = protein_path.stem

        # Proteinオブジェクトの作成
        protein = Protein.create(path=protein_path, id=protein_id)

        # 無効なオプション（負の値）- 必須パラメータを全て指定
        invalid_options = {
            "loop_min": -5,  # 負の値
            "helix_min": 30,
            "fragment_min": 5,
            "domain_min": 50,
            "pLDDT_min": 0,
            "local_contact_range": 5,
            "single_out": True,
            "domain_out": True,
        }

        # エラーが発生することを確認
        with pytest.raises(Exception):
            service.segment(protein, invalid_options, temp_output_dir)

    @pytest.mark.integration
    @pytest.mark.parametrize("protein_index", [0])
    def test_protein_segmentation_with_minimal_options(self, service, temp_output_dir, test_proteins, protein_index):
        """
        最小限のオプションでのテスト。

        Args:
            service: ProteinSegmentationServiceインスタンス
            temp_output_dir: 一時出力ディレクトリ
            test_proteins: テスト用のタンパク質ファイルのパスリスト
            protein_index: テスト対象のタンパク質のインデックス
        """
        # テスト対象のタンパク質ファイル
        protein_path = test_proteins[protein_index]
        protein_id = protein_path.stem

        # Proteinオブジェクトの作成
        protein = Protein.create(path=protein_path, id=protein_id)

        # AlphaCutterの必須パラメータを全て指定
        minimal_options = {
            "loop_min": 20,
            "helix_min": 30,
            "fragment_min": 5,
            "domain_min": 50,
            "pLDDT_min": 0,
            "local_contact_range": 5,
            "single_out": True,
            "domain_out": True,
        }

        # タンパク質セグメンテーションの実行
        segmented_proteins = service.segment(protein, minimal_options, temp_output_dir)

        # 結果の検証
        assert len(segmented_proteins) > 0, "セグメント化されたタンパク質が生成されていません"
