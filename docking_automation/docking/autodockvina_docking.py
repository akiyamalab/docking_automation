from __future__ import annotations

import tempfile
from pathlib import Path
from typing import TYPE_CHECKING, Any, Dict, List, Optional, Protocol, Tuple, Union

if TYPE_CHECKING:
    from docking_automation.infrastructure.repositories.hdf5_docking_result_repository import HDF5DockingResultRepository

import numpy as np
import numpy.typing as npt


class VinaProtocol(Protocol):
    """Vinaクラスのプロトコル定義"""

    def __init__(self, cpu: int = 1, seed: int = 0, verbosity: int = 1) -> None: ...
    def set_receptor(self, receptor_path: str) -> None: ...
    def set_ligand_from_file(self, ligand_path: str) -> None: ...
    def compute_vina_maps(self, center: List[float], box_size: List[float]) -> None: ...
    def write_maps(self, map_prefix_filename: str = "receptor", overwrite: bool = False) -> None: ...
    def load_maps(self, map_prefix_filename: str) -> None: ...
    def dock(self, exhaustiveness: int = 8, n_poses: int = 9, min_rmsd: float = 1.0) -> None: ...
    def write_poses(self, output_path: str, n_poses: int = 9, overwrite: bool = False) -> None: ...
    def energies(self) -> npt.NDArray[np.float64]: ...


from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein

from ..converters.molecule_converter import MoleculeConverter
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from .screening_tool import ScreeningTool
from .docking_parameters import DockingParameters, SpecificDockingParametersABC
from .docking_result import DockingResult
from vina import Vina


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
class AutoDockVina(ScreeningTool):
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

    def dock(
        self,
        parameters: DockingParameters,
        repository: Optional[HDF5DockingResultRepository] = None,
        verbose: bool = False
    ) -> List[DockingResult]:
        """
        AutoDock Vinaを使ってドッキング計算を実施する。

        複数の化合物に対してドッキング計算を行い、結果のリストを返す。
        CompoundSetにインデックス範囲が設定されている場合は、その範囲内の化合物のみを処理する。
        リポジトリが指定されている場合は、可能な場合は既存の結果を再利用する。

        Args:
            parameters: ドッキングパラメータ
            repository: 結果を保存/取得するリポジトリ（指定しない場合は再利用しない）
            verbose: 詳細なログを出力するかどうか

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
        compound_set: PreprocessedCompoundSet = common_params.compound_set

        if protein.file_path is None:
            raise ValueError("PreprocessedProtein.file_path が未設定のため AutoDock Vina を実行できません")
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

                # リポジトリが指定されている場合、既存の結果を確認
                if repository is not None and repository._exists(
                    protein.content_hash,
                    compound_hash
                ):
                    # 既存の結果を取得
                    if verbose:
                        print(f"化合物 {idx+1}/{len(compounds_to_process)} の結果を再利用します")
                    
                    result = repository.load_by_hashes(
                        protein.content_hash,
                        compound_hash
                    )
                    
                    if result is not None:
                        results.append(result)
                        if verbose:
                            print(f"化合物 {idx+1}/{len(compounds_to_process)} の結果を再利用しました（スコア: {result.docking_score}）")
                        continue  # 次の化合物へ
                
                # 既存の結果がない場合、通常通りドッキング計算を実行
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

    # --- ScreeningTool (旧 CacheableReceptorDocking) 準拠 ---
    # Vina の `compute_vina_maps` は receptor + grid box に対し atom-type 別の
    # potential grid を計算する処理で、典型的に 5〜30 秒かかる。同一 receptor に
    # 多数の ligand をドッキングする N×M スクリーニングでは、この結果を一度だけ
    # 計算して `.map` ファイル群に書き出しておけば、以降の docking では
    # `load_maps` で高速再利用できる。`docking_automation.docking.screening_tool.ScreeningTool` 参照。

    def prepare_receptor_cache(
        self,
        protein: Protein,
        grid_box,
        out_cache: Path,
        force: bool = False,
        seed: int = 0,
    ) -> Path:
        """Vina の atom-type 別 map ファイル群を生成してキャッシュする。

        Args:
            protein: 受容体。まだ PDBQT 化されていなくても `_preprocess_protein` で変換される。
            grid_box: ドッキングボックス。
            out_cache: 出力パス prefix (拡張子なし)。生成ファイル:
                `{out_cache}.C.map`, `{out_cache}.N.map`, ... など atom type 数個。
            force: True なら既存 `.map` 群を上書き。
            seed: Vina 初期化用 seed。map 計算自体には影響しない。

        Returns:
            out_cache (= 生成した .map 群の prefix)。load_maps にそのまま渡せる。
        """
        out_cache = Path(out_cache)
        # Vina が生成する map ファイル名は atom type (C_H, N_A, O_D など) 依存で
        # 受容体によって組が変わる可能性があるため、特定名ではなく glob で確認。
        out_cache.parent.mkdir(parents=True, exist_ok=True)
        existing = list(out_cache.parent.glob(f'{out_cache.name}.*.map'))
        if existing and not force:
            return out_cache

        preprocessed = self._preprocess_protein(protein)
        if preprocessed.file_path is None:
            raise ValueError("PreprocessedProtein.file_path が未設定のため map を生成できません")

        v = Vina(cpu=1, seed=seed, verbosity=0)
        v.set_receptor(str(preprocessed.file_path))
        v.compute_vina_maps(
            center=[float(c) for c in grid_box.center],
            box_size=[float(s) for s in grid_box.size],
        )
        v.write_maps(map_prefix_filename=str(out_cache), overwrite=True)
        return out_cache

    def dock_with_cache(
        self,
        cache: Path,
        ligand_paths: List[Path],
        grid_box,
        protein_content_hash: str,
        compound_content_hashes: Optional[List[str]] = None,
        exhaustiveness: int = 8,
        num_modes: int = 1,
        seed: int = 0,
    ) -> List[DockingResult]:
        """キャッシュ済み map を再利用して多数 ligand を docking。

        Args:
            cache: `prepare_receptor_cache` が返した prefix (例: `/path/to/prefix`)。
                `{prefix}.C.map` など atom-type 別 map が存在する前提。
            ligand_paths: PDBQT 形式の ligand ファイルリスト。
            grid_box: docking box。map と一致していれば良い。
            protein_content_hash: 結果 `DockingResult.protein_content_hash`。
            compound_content_hashes: None の場合はファイル stem を使用。
            exhaustiveness, num_modes, seed: Vina パラメータ。

        Returns:
            各 ligand の best pose 1 件の DockingResult。
        """
        v = Vina(cpu=1, seed=seed, verbosity=0)
        v.load_maps(map_prefix_filename=str(cache))

        results: List[DockingResult] = []
        for idx, lig_path in enumerate(ligand_paths):
            ch = (
                compound_content_hashes[idx]
                if compound_content_hashes is not None
                else lig_path.stem
            )
            try:
                tmp_dir = Path(tempfile.mkdtemp())
                out_pdbqt = tmp_dir / f"{lig_path.stem}_out.pdbqt"
                out_sdf = tmp_dir / f"{lig_path.stem}_out.sdf"

                v.set_ligand_from_file(str(lig_path))
                v.dock(exhaustiveness=exhaustiveness, n_poses=num_modes, min_rmsd=1.0)
                v.write_poses(str(out_pdbqt), n_poses=num_modes, overwrite=True)
                self.converter.pdbqt_to_sdf(out_pdbqt, out_sdf)
                scores = v.energies()

                results.append(
                    DockingResult(
                        result_path=out_sdf,
                        protein_id='',
                        compound_set_id=lig_path.parent.name,
                        compound_index=idx,
                        docking_score=float(scores[0, 0]),
                        protein_content_hash=protein_content_hash,
                        compound_content_hash=ch,
                        compoundset_content_hash=ch,
                        metadata={'tool': 'AutoDock Vina', 'source': 'vina_cached_maps'},
                    )
                )
            except Exception as e:
                print(f"ligand {lig_path.name} の dock_with_cache 失敗: {e}")
                continue

        return results
