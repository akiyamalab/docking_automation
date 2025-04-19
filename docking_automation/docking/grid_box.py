from typing import Final, List, Optional, Tuple, Union

import numpy as np
import numpy.typing as npt
from rdkit import Chem

from docking_automation.converters.molecule_converter import MoleculeConverter
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein

# TODO: [DDD] 値オブジェクトとしての実装を強化する
# - dataclass(frozen=True)への変換を検討する
# - 不変条件のバリデーションを強化する（サイズが正の値であることなど）
# - __eq__と__hash__メソッドを実装して等価性比較を可能にする
# - 値オブジェクト同士の演算メソッド（例：2つのGridBoxの結合や交差判定）を追加する


# 値オブジェクト
class GridBox:
    """
    ドッキング計算の領域を表すオブジェクト

    TODO: [DDD] このクラスをdataclass(frozen=True)に変換し、不変性を確保する
    例:
    @dataclass(frozen=True)
    class GridBox:
        center: np.ndarray
        size: np.ndarray

        def __post_init__(self):
            # 不変条件の検証
            if len(self.center) != 3 or len(self.size) != 3:
                raise ValueError("center と size は3次元配列である必要があります")
            if any(s <= 0 for s in self.size):
                raise ValueError("size はすべて正の値である必要があります")
    """

    __center: npt.NDArray[np.float64]
    __size: npt.NDArray[np.float64]

    def __init__(
        self,
        center: Optional[Union[npt.NDArray[np.float64], Tuple[float, float, float]]] = None,
        center_x: Optional[float] = None,
        center_y: Optional[float] = None,
        center_z: Optional[float] = None,
        size: Optional[Union[npt.NDArray[np.float64], Tuple[float, float, float]]] = None,
        size_x: Optional[float] = None,
        size_y: Optional[float] = None,
        size_z: Optional[float] = None,
    ):
        """
        GridBoxオブジェクトを初期化する。

        Args:
            center: 中心座標（NDArrayまたは(x, y, z)のタプル）
            center_x: 中心のx座標
            center_y: 中心のy座標
            center_z: 中心のz座標
            size: サイズ（NDArrayまたは(x, y, z)のタプル）
            size_x: x方向のサイズ
            size_y: y方向のサイズ
            size_z: z方向のサイズ

        Examples:
            >>> # タプルを使った初期化
            >>> grid_box = GridBox(center=(10.0, 10.0, 10.0), size=(20.0, 20.0, 20.0))
            >>>
            >>> # 個別の座標を使った初期化
            >>> grid_box = GridBox(
            ...     center_x=10.0, center_y=10.0, center_z=10.0,
            ...     size_x=20.0, size_y=20.0, size_z=20.0
            ... )
        """
        # 中心座標の初期化
        if center is not None:
            if isinstance(center, tuple):
                self.__center = np.array(center, dtype=np.float64)
            else:
                self.__center = center.copy()
        elif center_x is not None and center_y is not None and center_z is not None:
            self.__center = np.array([center_x, center_y, center_z], dtype=np.float64)
        else:
            raise ValueError(
                "中心座標が指定されていません。centerまたはcenter_x, center_y, center_zを指定してください。"
            )

        # サイズの初期化
        if size is not None:
            if isinstance(size, tuple):
                self.__size = np.array(size, dtype=np.float64)
            else:
                self.__size = size.copy()
        elif size_x is not None and size_y is not None and size_z is not None:
            self.__size = np.array([size_x, size_y, size_z], dtype=np.float64)
        else:
            raise ValueError("サイズが指定されていません。sizeまたはsize_x, size_y, size_zを指定してください。")

    @property
    def center(self) -> npt.NDArray[np.float64]:
        """
        中心座標を取得する。

        Returns:
            中心座標
        """
        return self.__center.copy()

    @property
    def size(self) -> npt.NDArray[np.float64]:
        """
        サイズを取得する。

        Returns:
            サイズ
        """
        return self.__size.copy()

    def __str__(self) -> str:
        """
        GridBoxの文字列表現を返す。

        Returns:
            GridBoxの文字列表現
        """
        return f"GridBox: 中心座標=({self.format_center()}), サイズ=({self.format_size()})"

    def __repr__(self) -> str:
        """
        GridBoxの再現可能な文字列表現を返す。

        Returns:
            GridBoxの再現可能な文字列表現
        """
        return f"GridBox(center={self.center}, size={self.size})"

    def format_center(self) -> str:
        """
        中心座標をフォーマットした文字列を返す。

        Returns:
            フォーマットされた中心座標の文字列
        """
        return f"{self.__center[0]:.3f}, {self.__center[1]:.3f}, {self.__center[2]:.3f}"

    def format_size(self) -> str:
        """
        サイズをフォーマットした文字列を返す。

        Returns:
            フォーマットされたサイズの文字列
        """
        return f"{self.__size[0]:.3f}, {self.__size[1]:.3f}, {self.__size[2]:.3f}"

    @classmethod
    def from_fpocket(cls, protein: Protein, pocket_rank: int = 1, keep_temp_files: bool = False) -> "GridBox":
        """
        fpocketを使用してタンパク質のポケット位置を予測し、GridBoxを生成する。

        Args:
            protein: 対象のタンパク質
            pocket_rank: ポケットのランク（1が最も有望）
            keep_temp_files: 一時ファイルを保持するかどうか。デバッグ目的で使用。

        Returns:
            予測されたGridBoxオブジェクト

        Raises:
            ValueError: fpocketの実行に失敗した場合や、指定されたランクのポケットが見つからない場合
        """
        # 循環インポートを避けるため、実行時にインポート
        from docking_automation.docking.fpocket_grid_box_predictor import (
            FpocketGridBoxPredictor,
        )

        # FpocketGridBoxPredictorを内部実装として使用
        predictor = FpocketGridBoxPredictor(keep_temp_files=keep_temp_files)
        return predictor.predict(protein, pocket_rank)

    @classmethod
    def from_fpocket_all(cls, protein: Protein, max_pockets: int = 3, keep_temp_files: bool = False) -> List["GridBox"]:
        """
        fpocketを使用してタンパク質の複数のポケット位置を予測し、GridBoxのリストを生成する。

        Args:
            protein: 対象のタンパク質
            max_pockets: 予測する最大ポケット数
            keep_temp_files: 一時ファイルを保持するかどうか。デバッグ目的で使用。

        Returns:
            予測されたGridBoxオブジェクトのリスト

        Raises:
            ValueError: fpocketの実行に失敗した場合
        """
        # 循環インポートを避けるため、実行時にインポート
        from docking_automation.docking.fpocket_grid_box_predictor import (
            FpocketGridBoxPredictor,
        )

        # FpocketGridBoxPredictorを内部実装として使用
        predictor = FpocketGridBoxPredictor(keep_temp_files=keep_temp_files)
        return predictor.predict_all(protein, max_pockets)

    @classmethod
    def from_crystal_ligand(cls, crystal_ligand: CompoundSet) -> "GridBox":
        """
        結晶構造リガンドからGridBoxを生成する。
        CompoundSetの1つめの構造を使用する。

        Args:
            crystal_ligand: 結晶構造リガンドの構造データ

        Returns:
            GridBoxオブジェクト
        """
        raise NotImplementedError("from_crystal_ligandメソッドは未実装です。")

    @classmethod
    def _eboxsize(cls, mol: Chem.Mol) -> int:
        """
        eBoxSizeアルゴリズムを使用してボックスサイズを計算する。
        https://www.brylinski.org/eboxsize

        Args:
            mol: RDKitの分子オブジェクト

        Returns:
            ボックスサイズ（偶数）
        """
        _GY_BOX_RATIO: Final[float] = 0.23

        # 重原子（水素以外の原子）の座標を取得
        coords = []
        conf = mol.GetConformer()
        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() > 1:  # 水素以外の原子
                # 型チェックを無視
                pos = conf.GetAtomPosition(int(atom_idx))  # type: ignore
                coords.append([pos.x, pos.y, pos.z])

        # NumPy配列に変換
        coords_array = np.array(coords)

        # 中心座標を計算
        center = coords_array.mean(axis=0)

        # 回転半径の二乗を計算
        sq_gyration = ((coords_array - center) ** 2).sum(axis=1).mean()

        # ボックスサイズを計算
        size = sq_gyration**0.5 / _GY_BOX_RATIO

        # 偶数に丸める（切り上げ）
        return int((size + 1) / 2) * 2

    @classmethod
    def from_compound(cls, compound: CompoundSet) -> "GridBox":
        """
        化合物からGridBoxを生成する。
        CompoundSetの1つめの構造を使用する。

        Args:
            compound: 化合物の構造データ

        Returns:
            GridBoxオブジェクト
        """
        # MoleculeConverterのインスタンスを作成
        converter = MoleculeConverter()

        # CompoundSetからRDKitの分子オブジェクトに変換
        mols = converter.compound_to_rdkit(compound)

        # 最初の分子を使用
        if not mols:
            raise ValueError("化合物が見つかりません")
        mol = mols[0]

        # 重原子（水素以外の原子）の座標を取得
        coords = []
        conf = mol.GetConformer()
        for atom_idx in range(mol.GetNumAtoms()):
            atom = mol.GetAtomWithIdx(atom_idx)
            if atom.GetAtomicNum() > 1:  # 水素以外の原子
                # 型チェックを無視
                pos = conf.GetAtomPosition(int(atom_idx))  # type: ignore
                coords.append([pos.x, pos.y, pos.z])

        # NumPy配列に変換
        coords_array = np.array(coords)

        # 中心座標を計算
        center = coords_array.mean(axis=0)

        # eBoxSizeアルゴリズムでボックスサイズを計算
        box_size = cls._eboxsize(mol)

        # すべての方向に同じサイズを使用
        size = np.array([box_size, box_size, box_size], dtype=np.float64)

        # GridBoxオブジェクトを作成して返す
        return cls(center=center, size=size)
