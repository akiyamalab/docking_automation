import glob
import importlib.resources
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

from docking_automation.molecule.protein import Protein
from docking_automation.molecule.protein_events import ProteinSegmented


class ProteinSegmentationService:
    """
    タンパク質のセグメンテーション（ドメイン分割・不要残基削除）を行うサービス。

    AlphaCutterスクリプトを使用して、タンパク質の低信頼度領域や非局所的相互作用を持たない
    フレキシブルな領域を除去し、構造的に安定したドメインを抽出・分割する。
    """

    def __init__(self, progress_callback: Optional[Callable[[str], None]] = None):
        """
        ProteinSegmentationServiceを初期化する。

        Args:
            progress_callback: 進捗状況を報告するコールバック関数（オプション）
        """
        self._progress_callback = progress_callback

    def _report_progress(self, message: str) -> None:
        """
        進捗状況を報告する。

        Args:
            message: 進捗メッセージ
        """
        if self._progress_callback:
            self._progress_callback(message)
        else:
            print(f"[ProteinSegmentation] {message}")

    def segment(self, protein: Protein, options: Dict[str, Any], final_output_dir: Path) -> List[Protein]:
        """
        タンパク質のセグメンテーションを実行する。

        Args:
            protein: 入力となる `Protein` オブジェクト
            options: AlphaCutterに渡すコマンドラインオプション（辞書形式）
            final_output_dir: 最終的な出力ファイル（分割されたPDBなど）を保存するディレクトリパス

        Returns:
            生成された新しい `Protein` オブジェクトのリスト
        """
        # スクリプトパスの取得
        script_path = self._get_script_path()

        # 出力ディレクトリの作成
        os.makedirs(final_output_dir, exist_ok=True)

        self._report_progress(f"タンパク質 {protein.id} のセグメンテーションを開始します")

        # 一時ディレクトリで処理を実行
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_dir_path = Path(temp_dir)

            # 入力ファイルを一時ディレクトリにコピー
            input_file_path = temp_dir_path / protein.path.name
            shutil.copy(protein.path, input_file_path)
            self._report_progress(f"入力ファイルを一時ディレクトリにコピーしました: {input_file_path}")

            # コマンドの構築
            command = self._build_alpha_cutter_command(script_path, input_file_path, options)
            self._report_progress("AlphaCutterコマンドを構築しました")

            # AlphaCutterの実行
            self._report_progress("AlphaCutterを実行中...")
            self._run_alpha_cutter(command, temp_dir_path)
            self._report_progress("AlphaCutterの実行が完了しました")

            # 出力ファイルの収集
            self._report_progress("出力ファイルを収集中...")
            output_file_paths = self._collect_outputs(temp_dir_path, final_output_dir)
            self._report_progress(f"{len(output_file_paths)}個の出力ファイルを収集しました")

            # 新しいProteinオブジェクトの生成
            self._report_progress("新しいProteinオブジェクトを生成中...")
            segmented_proteins = self._create_protein_objects(output_file_paths, protein)
            self._report_progress(f"{len(segmented_proteins)}個のセグメントが生成されました")

            # ドメインイベントの発行
            if segmented_proteins:
                segmented_protein_ids = [p.id for p in segmented_proteins]
                event = ProteinSegmented(original_protein_id=protein.id, segmented_protein_ids=segmented_protein_ids)
                # 最初のタンパク質にイベントを登録（他のタンパク質は既にProteinRegisteredイベントを持っている）
                segmented_proteins[0]._register_domain_event(event)
                self._report_progress("ドメインイベントを発行しました")

            return segmented_proteins

    def _get_script_path(self) -> Path:
        """
        AlphaCutterスクリプトの絶対パスを取得する。

        Returns:
            スクリプトの絶対パス
        """
        # importlib.resourcesを使用してパッケージ内のリソースパスを取得
        script_path = importlib.resources.files("docking_automation.infrastructure.external_tools").joinpath(
            "alpha_cutter.py"
        )
        return Path(str(script_path))

    def _build_alpha_cutter_command(
        self, script_path: Path, input_file_path: Path, options: Dict[str, Any]
    ) -> List[str]:
        """
        AlphaCutterコマンドを構築する。

        Args:
            script_path: AlphaCutterスクリプトのパス
            input_file_path: 入力ファイルのパス
            options: コマンドラインオプション（辞書形式）

        Returns:
            コマンドのリスト
        """
        command = ["python", str(script_path)]

        # 入力ファイルのパスを追加
        # AlphaCutterは現在のディレクトリのPDBファイルを処理するため、
        # 入力ファイルのパスは指定しませんが、カレントディレクトリを
        # 一時ディレクトリに設定することで、指定したファイルのみを処理します

        # オプションの追加
        for key, value in options.items():
            if isinstance(value, bool) and value:
                # フラグオプション（値なし）
                command.append(f"--{key}")
            elif not isinstance(value, bool):
                # 値を持つオプション
                command.append(f"--{key}")
                command.append(str(value))

        return command

    def _run_alpha_cutter(self, command: List[str], cwd: Path) -> None:
        """
        AlphaCutterを実行する。

        Args:
            command: 実行するコマンド
            cwd: カレントワーキングディレクトリ
        """
        try:
            result = subprocess.run(command, cwd=str(cwd), check=True, capture_output=True, text=True)
            # デバッグ用に標準出力と標準エラー出力を表示
            self._report_progress(f"AlphaCutter stdout: {result.stdout}")
            if result.stderr:
                self._report_progress(f"AlphaCutter stderr: {result.stderr}")
        except subprocess.CalledProcessError as e:
            self._report_progress(f"AlphaCutter実行エラー: {e}")
            self._report_progress(f"stdout: {e.stdout}")
            self._report_progress(f"stderr: {e.stderr}")
            raise

    def _collect_outputs(self, temp_dir: Path, final_output_dir: Path) -> List[Path]:
        """
        出力ファイルを収集し、最終出力ディレクトリに移動する。

        Args:
            temp_dir: 一時ディレクトリのパス
            final_output_dir: 最終出力ディレクトリのパス

        Returns:
            移動した出力ファイルのパスリスト
        """
        # 出力ファイルのパターン
        patterns = [
            "*_domain*.pdb",  # ドメイン分割されたPDBファイル
            "*_domains.pdb",  # ドメイン分割されたPDBファイル（別形式）
            "AFCT-OUT_*.csv",  # サマリーファイル
        ]

        output_file_paths = []

        for pattern in patterns:
            for file_path in temp_dir.glob(pattern):
                # 出力ディレクトリに移動
                dest_path = final_output_dir / file_path.name
                shutil.move(file_path, dest_path)

                # PDBファイルのみをリストに追加
                if dest_path.suffix.lower() == ".pdb":
                    output_file_paths.append(dest_path)

        return output_file_paths

    def _create_protein_objects(self, file_paths: List[Path], original_protein: Protein) -> List[Protein]:
        """
        ファイルパスのリストから新しいProteinオブジェクトを作成する。

        Args:
            file_paths: PDBファイルのパスリスト
            original_protein: 元のProteinオブジェクト

        Returns:
            作成されたProteinオブジェクトのリスト
        """
        proteins = []

        for path in file_paths:
            # ファイル名をベースにIDを生成（元のタンパク質IDを接頭辞として使用）
            protein_id = f"{original_protein.id}_{path.stem}"

            # 新しいProteinオブジェクトを作成
            protein = Protein.create(path=path, id=protein_id)
            proteins.append(protein)

        return proteins

    def get_segmentation_summary(self, segmented_proteins: List[Protein], original_protein: Protein) -> str:
        """
        セグメンテーション結果のサマリーを文字列として返す。

        Args:
            segmented_proteins: セグメンテーションされたタンパク質のリスト
            original_protein: 元のタンパク質

        Returns:
            セグメンテーション結果のサマリー文字列
        """
        if not segmented_proteins:
            return f"タンパク質 {original_protein.id} のセグメンテーション結果: セグメントなし"

        summary = [
            f"タンパク質 {original_protein.id} のセグメンテーション結果:",
            f"元のタンパク質: {original_protein}",
            f"セグメント数: {len(segmented_proteins)}個",
            "セグメント一覧:",
        ]

        for i, seg_protein in enumerate(segmented_proteins, 1):
            summary.append(f"  セグメント {i}: {seg_protein}")

        return "\n".join(summary)
