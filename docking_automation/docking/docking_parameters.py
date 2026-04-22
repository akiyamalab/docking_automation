from abc import ABC
from dataclasses import dataclass, field

from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein


# 値オブジェクト
class CommonDockingParameters:
    """
    全てのドッキングツールにおいて共通しているパラメータを保持するクラス。
    """

    def __init__(
        self,
        protein: PreprocessedProtein,
        compound_set: PreprocessedCompoundSet,
        grid_box,
    ):
        """
        CommonDockingParametersオブジェクトを初期化する。

        Args:
            protein: 前処理済みのタンパク質
            compound_set: 前処理済みの化合物セット
            grid_box: ドッキング計算の領域
        """
        self.protein = protein
        self.compound_set = compound_set
        self.grid_box = grid_box


# 値オブジェクト
class SpecificDockingParametersABC(ABC):
    """
    各ツール固有のパラメータを保持するクラス。
    """

    ...


@dataclass
class UniDockParameters(SpecificDockingParametersABC):
    """Uni-Dock GPU バッチドッキング固有パラメータ。"""
    search_mode: str = "balance"  # fast / balance / detail
    scoring: str = "vina"         # vina / vinardo
    num_modes: int = 1
    seed: int = 1
    verbosity: int = 0
    score_threshold_max: float = 5.0   # これ以上はペナルティスコアと判断
    score_threshold_min: float = -30.0  # これ以下も異常値と判断


# 値オブジェクト
class DockingParameters:
    """
    ドッキングツールに関するパラメータを保持するクラス。
    汎用的な部分と、ツール固有の部分を分けて保持する。
    """

    def __init__(self, common: CommonDockingParameters, specific: SpecificDockingParametersABC):
        self.common = common
        self.specific = specific
