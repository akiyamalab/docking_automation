import tempfile
from pathlib import Path
from typing import Any, Dict, List, Optional, Protocol, Tuple, Union

import numpy as np
import numpy.typing as npt


class VinaProtocol(Protocol):
    """Vinaクラスのプロトコル定義"""

    def __init__(self, cpu: int = 1, seed: int = 0, verbosity: int = 1) -> None: ...
    def set_receptor(self, receptor_path: str) -> None: ...
    def set_ligand_from_file(self, ligand_path: str) -> None: ...
    def compute_vina_maps(self, center: List[float], box_size: List[float]) -> None: ...
    def dock(self, exhaustiveness: int = 8, n_poses: int = 9, min_rmsd: float = 1.0) -> None: ...
    def write_poses(self, output_path: str, n_poses: int = 9, overwrite: bool = False) -> None: ...
    def energies(self) -> npt.NDArray[np.float64]: ...


from vina import Vina

from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein

from ..converters.molecule_converter import MoleculeConverter
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from .docking import DockingToolABC
from .docking_parameters import DockingParameters, SpecificDockingParametersABC
from .docking_result import DockingResult


# 値オブジェクト
class AutoDockVinaParameters(SpecificDockingParametersABC):
    """
    AutoDock Vina 固有のパラメータを保持するクラス。
    """

    def __init__(
        self,
        exhaustiveness: int = 8,
        num_modes: int = 9,
        seed: Optional[int] = None,
        max_compounds: Optional[int] = None,
        **kwargs: Any,
    ) -> None:
        """
        AutoDockVinaParametersオブジェクトを初期化する。

        Args:
            exhaustiveness: 探索の徹底度（デフォルト: 8）
            num_modes: 出力するポーズの数（デフォルト: 9）
            seed: 乱数シード（デフォルト: None）
            max_compounds: 処理する化合物の最大数（デフォルト: None、すべての化合物を処理）
            **kwargs: その他のパラメータ
        """
        self.params = {
            "exhaustiveness": exhaustiveness,
            "num_modes": num_modes,
        }

        if seed is not None:
            self.params["seed"] = seed

        if max_compounds is not None:
            self.params["max_compounds"] = max_compounds

        # その他のパラメータを追加
        self.params.update(kwargs)


# インフラ
class AutoDockVina(DockingToolABC):
    """
    AutoDock Vina を使ったドッキング計算を行うクラス。
    """

    def __init__(self) -> None:
        """
        AutoDockVinaオブジェクトを初期化する。
        """
        # TODO: converter が外から見える必要はないはず。
        self.converter = MoleculeConverter()

    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        """
        タンパク質について、AutoDock Vina用の前処理を行う。

        Args:
            protein: 前処理するタンパク質

        Returns:
            前処理済みのタンパク質
        """
        # 一時ディレクトリを作成
        temp_dir = Path(tempfile.mkdtemp())

        # PDBQTファイルに変換
        pdbqt_path = temp_dir / f"{protein.id}.pdbqt"
        self.converter.protein_to_pdbqt(protein, pdbqt_path)

        # 前処理済みのタンパク質を返す
        return PreprocessedProtein(file_path=pdbqt_path)

    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        """
        化合物について、AutoDock Vina用の前処理を行う。

        Args:
            compound_set: 前処理する化合物セット

        Returns:
            前処理済みの化合物セット
        """
        # 一時ディレクトリを作成
        temp_dir = Path(tempfile.mkdtemp())

        # PDBQTファイルに変換（複数の化合物に対応）
        pdbqt_paths = self.converter.compound_to_pdbqt(compound_set, temp_dir)

        # 化合物のハッシュ値キャッシュを取得
        # プライベート属性にアクセスするためのハック
        # 通常はこのような方法は避けるべきだが、、、、
        compound_hash_cache = getattr(compound_set, "_CompoundSet__compound_hash_cache").copy()

        # 前処理済みの化合物セットを返す
        return PreprocessedCompoundSet(file_paths=pdbqt_paths, compound_hash_cache=compound_hash_cache)

    def dock(self, parameters: DockingParameters, verbose: bool = False) -> List[DockingResult]:
        """
        AutoDock Vinaを使ってドッキング計算を実施する。

        複数の化合物に対してドッキング計算を行い、結果のリストを返す。
        CompoundSetにインデックス範囲が設定されている場合は、その範囲内の化合物のみを処理する。

        Args:
            parameters: ドッキングパラメータ

        Returns:
            ドッキング結果のリスト
        """
        # パラメータを取得
        common_params = parameters.common
        specific_params = parameters.specific

        if not isinstance(specific_params, AutoDockVinaParameters):
            raise ValueError("specific_paramsはAutoDockVinaParametersのインスタンスである必要があります。")

        # 前処理済みのタンパク質と化合物セット
        protein = common_params.protein
        compound_set: Union[CompoundSet, PreprocessedCompoundSet] = common_params.compound_set
        grid_box = common_params.grid_box

        # ファイルパスを取得

        # 一時ディレクトリを作成
        temp_dir = Path(tempfile.mkdtemp())

        # グリッドボックスの中心とサイズを取得
        center = grid_box.center
        size = grid_box.size

        # 結果を格納するリスト
        results = []

        # 処理する化合物の最大数を取得
        max_compounds = specific_params.params.get("max_compounds")

        # インデックス範囲の初期化
        start_index = 0

        # CompoundSetのプロパティを取得して、インデックスリストまたはインデックス範囲が設定されているかどうかを確認
        try:
            properties = compound_set.get_properties()
            
            # インデックスリストが設定されている場合
            indices = properties.get("indices")
            if indices is not None:
                # インデックスリストが設定されている場合は、start_indexは0のままでOK
                # 実際のインデックスはcompound_indexの計算時に使用する
                pass
            # インデックス範囲が設定されている場合
            elif "index_range" in properties:
                index_range = properties["index_range"]
                start_index = index_range["start"]
        except Exception as e:
            print(f"インデックス情報の取得中にエラーが発生しました: {e}")

        # 化合物の数を取得（各タスクで処理する化合物数）
        task_compounds = len(compound_set.file_paths)

        # 処理する化合物の数を決定
        if max_compounds is not None and max_compounds > 0 and max_compounds < task_compounds:
            compounds_to_process = compound_set.file_paths[:max_compounds]
            if verbose:
                print(f"ドッキング計算を開始します（全{task_compounds}化合物中、最初の{max_compounds}化合物）...")
        else:
            compounds_to_process = compound_set.file_paths
            if verbose:
                print(f"ドッキング計算を開始します（全{task_compounds}化合物）...")

        # 各化合物に対してドッキング計算を実行
        for idx, compound_path in enumerate(compounds_to_process):
            try:
                if verbose:
                    print(f"化合物 {idx+1}/{len(compounds_to_process)} を処理中...")

                # 結果ファイルのパス（化合物ごとに区別）
                output_pdbqt = temp_dir / f"output_{idx}.pdbqt"
                output_sdf = temp_dir / f"output_{idx}.sdf"

                # Vinaオブジェクトを作成
                # 並列計算は外側でやるので内部は1スレッドで実行
                v = Vina(cpu=1, seed=1, verbosity=0)

                # 受容体を設定
                v.set_receptor(str(protein.file_path))  # タンパク質は1つのみ

                # リガンドを設定
                v.set_ligand_from_file(str(compound_path))

                # スコア関数を設定（グリッドボックスの中心と大きさを指定）
                v.compute_vina_maps(center=[center[0], center[1], center[2]], box_size=[size[0], size[1], size[2]])

                # ドッキング計算を実行
                v.dock(
                    exhaustiveness=specific_params.params.get("exhaustiveness", 8),
                    n_poses=specific_params.params.get("num_modes", 9),
                    min_rmsd=1.0,
                )

                # 結果を保存
                v.write_poses(str(output_pdbqt), n_poses=specific_params.params.get("num_modes", 9), overwrite=True)

                # 結果をSDFに変換
                self.converter.pdbqt_to_sdf(output_pdbqt, output_sdf)

                # スコアを取得
                scores: npt.NDArray[np.float64] = v.energies()

                # メタデータを作成
                metadata = {
                    "tool": "AutoDock Vina",
                    "parameters": specific_params.params,
                    "scores": scores,
                    "pose_path": str(output_sdf),
                }

                # 化合物のハッシュ値を取得
                compound_hash = compound_set.get_compound_hash(idx)

                # 実際の化合物インデックスを計算
                # インデックスリストが設定されている場合は、そのリスト内のインデックスを使用
                compound_index = idx
                if "indices" in properties:
                    # インデックスリストが設定されている場合は、そのリスト内のインデックスを使用
                    compound_index = properties["indices"][idx]
                else:
                    # インデックス範囲が設定されている場合は、start_indexを加算
                    compound_index = start_index + idx

                # DockingResultオブジェクトを作成
                result = DockingResult(
                    result_path=output_sdf,  # SDFファイルのパスを設定
                    protein_id=protein.file_path.stem,  # タンパク質は1つのみ
                    compound_set_id=compound_path.stem.split("_")[0],  # 化合物セットID（ファイル名から抽出）
                    compound_index=compound_index,  # 実際の化合物インデックス
                    docking_score=scores[0, 0],
                    protein_content_hash=protein.content_hash,
                    compound_content_hash=compound_hash,
                    compoundset_content_hash=compound_set.content_hash,
                    metadata=metadata,
                )

                results.append(result)
                if verbose:
                    print(
                        f"化合物 {idx+1}/{len(compounds_to_process)} のドッキング計算が完了しました（スコア: {scores[0,0]}）"
                    )

            except Exception as e:
                print(f"化合物 {idx+1}/{len(compounds_to_process)} の処理中にエラーが発生しました: {str(e)}")
                # エラーが発生しても処理を継続
                continue

        if verbose:
            print(f"ドッキング計算が完了しました（成功: {len(results)}/{len(compounds_to_process)}）")

        if not results:
            raise ValueError("有効なドッキング結果が得られませんでした。")

        return results
