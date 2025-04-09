import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any, List

# 外部ライブラリのインポート
import openbabel.openbabel as ob  # type: ignore[import-untyped]
import openbabel.pybel as pybel  # type: ignore[import-untyped]
from rdkit import Chem

from docking_automation.infrastructure.utilities.error_utils import capture_stderr
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein

# OpenBabelのwarningを抑制
ob.obErrorLog.SetOutputLevel(ob.obError)

# meekoをインポート
from meeko import (  # type: ignore[import-untyped]
    MoleculePreparation,
    PDBQTMolecule,
    PDBQTWriterLegacy,
    RDKitMolCreate,
)


class MoleculeConverter:
    """
    分子変換処理を行うクラス。

    OpenBabel、RDKit、meeko等を使用して、様々な分子形式間の変換を行う。
    """

    def protein_to_openbabel(self, protein: Protein, quiet: bool = False) -> pybel.Molecule:
        """
        ProteinオブジェクトからOpenBabelの分子オブジェクトに変換する。

        Args:
            protein: 変換対象のProteinオブジェクト
            quiet: Trueの場合、warningを抑制する

        Returns:
            OpenBabelの分子オブジェクト
        """
        # ファイル形式を判断
        # TODO: これは Protein がvalidationすべきではないか？
        file_format = protein.path.suffix.lstrip(".")
        if file_format.lower() != "pdb":
            raise ValueError(f"サポートされていないファイル形式です: {file_format}")

        try:
            if quiet:
                # 標準エラー出力をキャプチャ
                with capture_stderr():
                    # PDBファイルを読み込む
                    mol: pybel.Molecule = next(pybel.readfile(file_format, str(protein.path)))
                print(f"タンパク質ファイル {protein.path.name} を読み込みました。")
            else:
                # 標準エラー出力をそのまま表示
                mol = next(pybel.readfile(file_format, str(protein.path)))
            return mol
        except Exception as e:
            raise ValueError(f"タンパク質ファイル {protein.path.name} の読み込み中にエラーが発生しました")

    # TODO: id は与えられなくても良い。その場合は UUID のようなものが生成されるように Protein 側でなっている。
    def openbabel_to_protein(self, ob_mol: Any, id: str, output_path: Path) -> Protein:
        """
        OpenBabelの分子オブジェクトからProteinオブジェクトに変換する。

        Args:
            ob_mol: 変換対象のOpenBabel分子オブジェクト
            id: 生成するProteinのID
            output_path: 出力先のパス

        Returns:
            生成されたProteinオブジェクト
        """
        # 出力ディレクトリが存在しない場合は作成
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # PDBファイルとして保存
        ob_mol.write("pdb", str(output_path), overwrite=True)

        # Proteinオブジェクトを作成
        return Protein(output_path, id)

    def compound_to_rdkit(self, compound: CompoundSet) -> List[Chem.Mol]:
        """
        CompoundSetオブジェクトからRDKitの分子オブジェクトのリストに変換する。

        Args:
            compound: 変換対象のCompoundSetオブジェクト

        Returns:
            RDKitの分子オブジェクトのリスト
        """
        # ファイル形式を判断
        file_format = compound.path.suffix.lstrip(".")

        # gzipで圧縮されている場合は一時ファイルに展開
        # TODO: gzip で圧縮されているか否かをユーザが意識しないようにする 汎用 関数を作るべき。
        if compound.path.suffix.lower() == ".gz":
            # MEMO: わざわざ delete=False にする理由は何か？
            with tempfile.NamedTemporaryFile(delete=False, suffix=compound.path.stem) as temp_file:
                temp_path = Path(temp_file.name)
                with gzip.open(compound.path, "rb") as f_in:
                    with open(temp_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                # 実際のファイル形式を取得
                file_format = compound.path.stem.split(".")[-1].lower()

                if file_format == "sdf":
                    # SDFファイルを読み込む
                    suppl = Chem.SDMolSupplier(str(temp_path))
                    # TODO: 何件の化合物がファイルに含まれており、何件を正常に読み込めたかをログに残すべきである。
                    mols: list[Chem.Mol] = [mol for mol in suppl if mol is not None]
                else:
                    raise ValueError(f"サポートされていないファイル形式です: {file_format}")

                # 一時ファイルを削除
                os.unlink(temp_path)
        else:
            if file_format.lower() == "sdf":
                # SDFファイルを読み込む
                suppl = Chem.SDMolSupplier(str(compound.path))
                mols = [mol for mol in suppl if mol is not None]
            else:
                raise ValueError(f"サポートされていないファイル形式です: {file_format}")

        return mols

    # TODO: id は optional であるべき。
    def rdkit_to_compound(self, mols: List[Chem.Mol], id: str, output_path: Path) -> CompoundSet:
        """
        RDKitの分子オブジェクトのリストからCompoundSetオブジェクトに変換する。

        Args:
            mols: 変換対象のRDKit分子オブジェクトのリスト
            id: 生成するCompoundSetのID
            output_path: 出力先のパス

        Returns:
            生成されたCompoundSetオブジェクト
        """
        # 出力ディレクトリが存在しない場合は作成
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # SDFファイルとして保存
        writer = Chem.SDWriter(str(output_path))
        for mol in mols:
            writer.write(mol)
        writer.close()

        # CompoundSetオブジェクトを作成
        return CompoundSet(output_path, id)

    def protein_to_pdbqt(self, protein: Protein, output_path: Path) -> Path:
        """
        Proteinオブジェクトからpdbqtファイルに変換する。
        prepare_receptorコマンドを使用して変換を行う。

        Args:
            protein: 変換対象のProteinオブジェクト
            output_path: 出力先のパス

        Returns:
            生成されたpdbqtファイルのパス

        Raises:
            ValueError: prepare_receptorコマンドの実行に失敗した場合
        """
        # 出力ディレクトリが存在しない場合は作成
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # 一時ディレクトリを作成
        temp_dir = Path(tempfile.mkdtemp())

        try:
            # 水素原子を削除した一時ファイルを作成
            # ファイル形式を判断
            file_format = protein.path.suffix.lstrip(".")

            # pybelを使用して分子を読み込む
            try:
                # 標準エラー出力をキャプチャ
                with capture_stderr():
                    # 分子を読み込む
                    mol = next(pybel.readfile(file_format, str(protein.path)))

                    # 水素原子を削除
                    mol.OBMol.DeleteHydrogens()
            except Exception as e:
                raise ValueError(f"タンパク質ファイル {protein.path.name} の処理中にエラーが発生しました: {e}")

            # 水素原子を削除した一時ファイルのパス
            temp_protein_path = temp_dir / f"{protein.id}_no_hydrogens.pdb"

            # 水素原子を削除した分子をPDBファイルとして保存
            mol.write("pdb", str(temp_protein_path), overwrite=True)

            # prepare_receptorコマンドを使用してPDBQTに変換
            # prepare_receptorコマンドのオプション
            # -r: 入力ファイル
            # -o: 出力ファイル
            # -A: 修復タイプ（hydrogens: 水素を追加）
            # -U: クリーンアップタイプ（nphs: 非極性水素を削除）
            cmd = [
                "prepare_receptor",
                "-r",
                str(temp_protein_path),
                "-o",
                str(output_path),
                "-A",
                "hydrogens",
                "-U",
                "nphs",
            ]
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            return output_path
        except subprocess.CalledProcessError as e:
            # TODO: ログ出力は標準ライブラリを利用したい
            error_message = f"prepare_receptorの実行中にエラーが発生しました: {e}\n"
            error_message += f"標準エラー出力: {e.stderr}\n"
            error_message += f"標準出力: {e.stdout}\n"
            error_message += f"コマンド: {' '.join(cmd)}\n"
            error_message += f"終了コード: {e.returncode}"
            raise ValueError(error_message)
        except Exception as e:
            raise ValueError(f"タンパク質の前処理中にエラーが発生しました: {e}")
        finally:
            # 一時ディレクトリを削除
            shutil.rmtree(temp_dir)

    def _convert_to_pdbqt_with_meeko(self, mol: Chem.Mol, output_path: Path) -> Path:
        """
        RDKitの分子オブジェクトをmeekoを使用してPDBQTに変換する。

        Args:
            mol: 変換対象のRDKit分子オブジェクト
            output_path: 出力先のパス

        Returns:
            生成されたPDBQTファイルのパス
        """
        # 明示的な水素原子を追加
        mol_with_h = Chem.AddHs(mol)

        # meekoを使用してPDBQTに変換
        preparator = MoleculePreparation()
        # prepare()メソッドの戻り値を使用
        molecule_setups = preparator.prepare(mol_with_h)

        # molecule_setupsが空でないことを確認
        if not molecule_setups:
            raise ValueError("Meekoによる分子の準備に失敗しました。")

        # 最初のMoleculeSetupインスタンスを使用してPDBQTを生成
        pdbqt_string, is_ok, err_msg = PDBQTWriterLegacy.write_string(molecule_setups[0])
        if not is_ok:
            raise ValueError(f"PDBQTの生成に失敗しました: {err_msg}")

        # PDBQTファイルに保存
        with open(output_path, "w") as f:
            f.write(pdbqt_string)

        return output_path

    def compound_to_pdbqt(self, compound: CompoundSet, output_dir: Path, verbose: bool = False) -> List[Path]:
        """
        CompoundSetオブジェクトからpdbqtファイルに変換する。
        meekoを使用して変換を行う。

        複数の化合物を含むSDFファイルに対応し、各化合物に対して個別のPDBQTファイルを生成する。
        CompoundSetにインデックス範囲が設定されている場合は、その範囲内の化合物のみを処理する。

        Args:
            compound: 変換対象のCompoundSetオブジェクト
            output_dir: 出力先のディレクトリ

        Returns:
            生成されたpdbqtファイルのパスのリスト
        """
        # 出力ディレクトリが存在しない場合は作成
        output_dir.mkdir(parents=True, exist_ok=True)

        # 一時的なSDFファイルを作成
        with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as temp_file:
            temp_path = Path(temp_file.name)

            # ファイル形式を判断
            file_format = compound.path.suffix.lower()

            # gzipで圧縮されている場合は一時ファイルに展開
            if file_format == ".gz":
                with gzip.open(compound.path, "rb") as f_in:
                    with open(temp_path, "wb") as f_out:
                        shutil.copyfileobj(f_in, f_out)

                # 実際のファイル形式を取得
                actual_format = compound.path.stem.split(".")[-1].lower()
                input_path = temp_path
            else:
                input_path = compound.path

            try:
                # RDKitを使用して分子を読み込む
                if input_path.suffix.lower() == ".sdf":
                    # SDFファイルを読み込む
                    suppl = Chem.SDMolSupplier(str(input_path))

                    # 各化合物に対してPDBQTファイルを生成
                    pdbqt_paths = []
                    valid_count = 0
                    total_count = 0

                    # CompoundSetのプロパティを取得して、インデックス範囲が設定されているかどうかを確認
                    properties = compound.get_properties()
                    index_range = properties.get("index_range")

                    # インデックス範囲が設定されている場合は、その範囲内の化合物数を表示
                    if index_range is not None:
                        start_index = index_range["start"]
                        end_index = index_range["end"]
                        task_compounds = end_index - start_index
                        if verbose:
                            print(
                                f"化合物の前処理を開始します（インデックス範囲: {start_index}-{end_index-1}、化合物数: {task_compounds}）..."
                            )
                    else:
                        if verbose:
                            # インデックス範囲が設定されていない場合は、全体の化合物数を表示
                            print(f"化合物の前処理を開始します...")

                    # インデックス範囲が設定されている場合は、その範囲内の化合物のみを処理
                    for idx, mol in enumerate(suppl):
                        # インデックス範囲が設定されている場合は、その範囲内かどうかをチェック
                        if index_range is not None:
                            start_index = index_range["start"]
                            end_index = index_range["end"]
                            if idx < start_index or idx >= end_index:
                                continue

                        # 処理対象の化合物をカウント
                        total_count += 1

                        if mol is None:
                            print(f"警告: インデックス {idx} の化合物は無効です。スキップします。")
                            continue

                        try:
                            # 出力ファイルのパス
                            output_path = output_dir / f"{compound.id}_{idx}.pdbqt"

                            # meekoを使用してPDBQTに変換
                            self._convert_to_pdbqt_with_meeko(mol, output_path)

                            # ファイルが空の場合はエラー
                            if output_path.stat().st_size == 0:
                                print(f"警告: インデックス {idx} の化合物のPDBQTファイルが空です。スキップします。")
                                continue

                            pdbqt_paths.append(output_path)
                            valid_count += 1

                            if verbose and valid_count % 10 == 0:
                                print(f"進捗: {valid_count}/{total_count} 化合物を処理しました")

                        except Exception as e:
                            print(f"警告: インデックス {idx} の化合物の処理中にエラーが発生しました: {e}")
                            continue

                    if verbose:
                        print(f"化合物の前処理が完了しました（成功: {valid_count}/{total_count}）")

                    if not pdbqt_paths:
                        raise ValueError("有効なPDBQTファイルが生成されませんでした。")

                    return pdbqt_paths
                else:
                    raise ValueError(f"サポートされていないファイル形式です: {input_path.suffix.lower()}")
            except Exception as e:
                raise ValueError(f"PDBQTファイルの生成中にエラーが発生しました: {e}")
            finally:
                # 一時ファイルを削除
                if file_format == ".gz":
                    os.unlink(temp_path)

    def pdbqt_to_sdf(self, pdbqt_path: Path, output_path: Path) -> Path:
        """
        pdbqtファイルからsdfファイルに変換する。
        meekoとRDKitを使用して変換を行う。

        Args:
            pdbqt_path: 変換対象のpdbqtファイルのパス
            output_path: 出力先のパス

        Returns:
            生成されたsdfファイルのパス

        Raises:
            ValueError: 変換に失敗した場合

        Note:
            pdbqtファイルは、meekoによって作成されたものを想定しており、
            それ以外のファイルを変換する場合にはエラーが発生する。
        """
        # 出力ディレクトリが存在しない場合は作成
        output_path.parent.mkdir(parents=True, exist_ok=True)

        mols: List[Any] = self.pdbqt_to_rdkit(pdbqt_path)

        # SDFファイルとして保存
        with Chem.SDWriter(str(output_path)) as writer:
            for mol in mols:
                writer.write(mol)

        return output_path

    def pdbqt_to_rdkit(self, pdbqt_path: Path) -> List[Chem.Mol]:
        """
        pdbqtファイルからRDKitの分子オブジェクトのリストに変換する。

        Args:
            pdbqt_path: 変換対象のpdbqtファイルのパス

        Returns:
            RDKitの分子オブジェクトのリスト

        Note:
            1つのpdbqtファイルには、同じ分子構造の異なるconformerが含まれることがある。
            そのため、出力は list[Chem.Mol] となる。
        """

        with open(pdbqt_path, "r") as f:
            pdbqt_string = f.read()
        pdbqt_mol = PDBQTMolecule(pdbqt_string)
        rdkit_mol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        return rdkit_mol_list
