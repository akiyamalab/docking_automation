from typing import Optional, Tuple, Union
import numpy as np
import numpy.typing as npt

from docking_automation.molecule.compound_set import CompoundSet

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
        size_z: Optional[float] = None
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
            raise ValueError("中心座標が指定されていません。centerまたはcenter_x, center_y, center_zを指定してください。")
        
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

    @classmethod
    def from_crystal_ligand(cls, crystal_ligand: CompoundSet) -> 'GridBox':
        """
        結晶構造リガンドからGridBoxを生成する。
        CompoundSetの1つめの構造を使用する。
        
        Args:
            crystal_ligand: 結晶構造リガンドの構造データ
        
        Returns:
            GridBoxオブジェクト
        """
        raise NotImplementedError("from_crystal_ligandメソッドは未実装です。")