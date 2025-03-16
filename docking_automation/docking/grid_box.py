from pathlib import Path
from typing import Optional
import numpy as np
import numpy.typing as npt


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
        ...
