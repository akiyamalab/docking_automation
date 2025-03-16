from typing import Optional
import numpy as np
import numpy.typing as npt

# 値オブジェクト
class GridBox:
    """
    ドッキング計算の領域を表すオブジェクト
    """
    __center: npt.NDArray[np.float64]
    __size: npt.NDArray[np.float64]
    
    def __init__(
        self,
        center: Optional[npt.NDArray[np.float64]] = None,
        center_x: Optional[float] = None,
        center_y: Optional[float] = None,
        center_z: Optional[float] = None,
        size: Optional[npt.NDArray[np.float64]] = None,
        size_x: Optional[float] = None,
        size_y: Optional[float] = None,
        size_z: Optional[float] = None
    ):
        """
        GridBoxオブジェクトを初期化する。
        
        Args:
            center: 中心座標
            center_x: 中心のx座標
            center_y: 中心のy座標
            center_z: 中心のz座標
            size: サイズ
            size_x: x方向のサイズ
            size_y: y方向のサイズ
            size_z: z方向のサイズ
        """
        raise NotImplementedError()
    
    def get_center(self) -> npt.NDArray[np.float64]:
        """
        中心座標を取得する。
        
        Returns:
            中心座標
        """
        raise NotImplementedError()
    
    def get_size(self) -> npt.NDArray[np.float64]:
        """
        サイズを取得する。
        
        Returns:
            サイズ
        """
        raise NotImplementedError()
