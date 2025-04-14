import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

# 循環インポートを避けるため、型ヒントではAnyを使用
from docking_automation.molecule.protein import Protein


class FpocketGridBoxPredictor:
    """
    fpocketを使用してタンパク質のポケット位置を予測し、GridBoxを生成するクラス。

    fpocketはタンパク質のPDBファイルを解析し、ポケット（リガンド結合部位）を検出するツールです。
    検出されたポケットの情報を基に、ドッキング計算に使用するGridBoxを生成します。
    """

    def __init__(self, keep_temp_files: bool = False):
        """
        FpocketGridBoxPredictorを初期化する。

        Args:
            keep_temp_files: 一時ファイルを保持するかどうか。デバッグ目的で使用。
        """
        self.__keep_temp_files = keep_temp_files

    def predict(self, protein: Protein, pocket_rank: int = 1) -> Any:
        """
        指定されたランクのポケットに基づいてGridBoxを予測する。

        Args:
            protein: 対象のタンパク質
            pocket_rank: ポケットのランク（1が最も有望）

        Returns:
            予測されたGridBox

        Raises:
            ValueError: fpocketの実行に失敗した場合や、指定されたランクのポケットが見つからない場合
        """
        # predict_allを使用して効率的に実装
        grid_boxes = self.predict_all(protein, max_pockets=pocket_rank)

        # 指定されたランクのポケットが存在するか確認
        if len(grid_boxes) < pocket_rank:
            raise ValueError(f"ランク {pocket_rank} のポケットが見つかりません。利用可能なランク数: {len(grid_boxes)}")

        # 指定されたランクのポケットのGridBoxを返す
        return grid_boxes[pocket_rank - 1]

    def predict_all(self, protein: Protein, max_pockets: int = 3) -> List[Any]:
        """
        複数のポケットに基づいてGridBoxのリストを予測する。

        Args:
            protein: 対象のタンパク質
            max_pockets: 予測する最大ポケット数

        Returns:
            予測されたGridBoxのリスト

        Raises:
            ValueError: fpocketの実行に失敗した場合
        """
        # 一時ディレクトリを作成
        with tempfile.TemporaryDirectory() as temp_dir_str:
            temp_dir = Path(temp_dir_str)

            try:
                # fpocketを実行
                output_pdb_path = self._run_fpocket(protein.path, temp_dir)

                # fpocketの出力を解析
                pocket_coords = self._parse_fpocket_output(output_pdb_path)

                # 利用可能なポケットのランクを取得
                available_ranks = sorted(pocket_coords.keys())

                # 指定された最大数まで、または利用可能なすべてのポケットについてGridBoxを作成
                grid_boxes = []
                for rank in available_ranks[:max_pockets]:
                    grid_box = self._create_grid_box_from_pocket(pocket_coords[rank])
                    grid_boxes.append(grid_box)

                return grid_boxes

            finally:
                # fpocketの出力ディレクトリのパスを取得
                output_name = protein.path.stem
                fpocket_output_dir = protein.path.parent / f"{output_name}_out"

                # 一時ファイルを保持する場合
                if self.__keep_temp_files:
                    # fpocketの出力ディレクトリをコピー
                    if fpocket_output_dir.exists():
                        target_dir = protein.path.parent / f"{output_name}_fpocket_debug"
                        shutil.copytree(fpocket_output_dir, target_dir, dirs_exist_ok=True)
                        print(f"デバッグ用にfpocketの出力を保存しました: {target_dir}")
                else:
                    # 一時ファイルを保持しない場合は、fpocketの出力ディレクトリを削除
                    if fpocket_output_dir.exists():
                        shutil.rmtree(fpocket_output_dir)

    def _run_fpocket(self, protein_path: Path, output_dir: Path) -> Path:
        """
        fpocketを実行し、結果のパスを返す。

        Args:
            protein_path: タンパク質のPDBファイルパス
            output_dir: 出力ディレクトリ

        Returns:
            fpocketの出力PDBファイルのパス

        Raises:
            ValueError: fpocketの実行に失敗した場合
        """
        # 出力ディレクトリが存在しない場合は作成
        output_dir.mkdir(parents=True, exist_ok=True)

        # fpocketコマンドを構築
        cmd = [
            "fpocket",
            "-f",
            str(protein_path.relative_to(Path.cwd())),  # 絶対パスを使用
        ]

        try:
            # fpocketを実行
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            # 出力ファイルのパスを構築
            # fpocketは入力ファイル名_outディレクトリに結果を出力する
            output_name = protein_path.stem
            fpocket_output_dir = protein_path.parent / f"{output_name}_out"
            output_pdb_path = fpocket_output_dir / f"{output_name}_out.pdb"

            if not output_pdb_path.exists():
                raise ValueError(f"fpocketの出力ファイル {output_pdb_path} が見つかりません")

            return output_pdb_path

        except subprocess.CalledProcessError as e:
            error_message = f"fpocketの実行中にエラーが発生しました: {e}\n"
            error_message += f"標準エラー出力: {e.stderr}\n"
            error_message += f"標準出力: {e.stdout}\n"
            error_message += f"コマンド: {' '.join(cmd)}\n"
            error_message += f"終了コード: {e.returncode}"
            raise ValueError(error_message)

    def _parse_fpocket_output(self, output_pdb_path: Path) -> Dict[int, List[Tuple[float, float, float]]]:
        """
        fpocketの出力PDBファイルを解析し、ポケットごとの座標データを返す。

        Args:
            output_pdb_path: fpocketの出力PDBファイルのパス

        Returns:
            ポケット番号をキー、座標のリストを値とする辞書
        """
        pocket_coords: Dict[int, List[Tuple[float, float, float]]] = {}

        # PDBファイルを読み込む
        with open(output_pdb_path, "r") as f:
            for line in f:
                # STP残基の行を検出
                if line.startswith("HETATM") and "STP" in line:
                    # PDBフォーマットから情報を抽出
                    # 例: HETATM    1  O   STP     1      40.267  27.857  40.665  0.00  0.00          O
                    try:
                        residue_number = int(line[22:26].strip())  # 残基番号
                        x = float(line[30:38].strip())  # X座標
                        y = float(line[38:46].strip())  # Y座標
                        z = float(line[46:54].strip())  # Z座標

                        # 残基番号ごとに座標を収集
                        if residue_number not in pocket_coords:
                            pocket_coords[residue_number] = []

                        pocket_coords[residue_number].append((x, y, z))

                    except (ValueError, IndexError) as e:
                        print(f"警告: PDB行の解析中にエラーが発生しました: {line.strip()}, エラー: {e}")
                        continue

        if not pocket_coords:
            raise ValueError(f"STP残基が見つかりません: {output_pdb_path}")

        return pocket_coords

    def _create_grid_box_from_pocket(self, pocket_coords: List[Tuple[float, float, float]]) -> Any:
        """
        ポケットの座標データからGridBoxを作成する。

        Args:
            pocket_coords: ポケットの座標データのリスト

        Returns:
            作成されたGridBox
        """
        # 座標をNumPy配列に変換
        coords_array = np.array(pocket_coords)

        # 中心座標を計算（平均）
        center = coords_array.mean(axis=0)

        # 回転半径の二乗を計算
        sq_gyration = ((coords_array - center) ** 2).sum(axis=1).mean()

        # ボックスサイズを計算（eboxsizeアルゴリズム）
        _GY_BOX_RATIO = 0.23
        size = sq_gyration**0.5 / _GY_BOX_RATIO

        # 偶数に丸める（切り上げ）
        box_size = int((size + 1) / 2) * 2

        # サイズが10未満の場合は10に調整
        if box_size < 10:
            print(f"警告: ドッキング範囲のサイズが小さすぎるため、最小サイズ10に調整します（元のサイズ: {box_size}）")
            box_size = 10

        # すべての方向に同じサイズを使用
        size = np.array([box_size, box_size, box_size], dtype=np.float64)

        # GridBoxオブジェクトを作成して返す
        # 循環インポートを避けるため、実行時にインポート
        from docking_automation.docking.grid_box import GridBox

        return GridBox(center=center, size=size)
