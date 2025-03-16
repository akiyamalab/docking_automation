from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from .docking_result import DockingResult
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from .docking_parameters import DockingParameters, SpecificDockingParametersABC
from .docking import DockingToolABC

# 値オブジェクト
class REstrettoParameters(SpecificDockingParametersABC):
    """
    REstretto 固有のパラメータを保持するクラス。
    """
    def __init__(self, **kwargs):
        """
        REstrettoParametersオブジェクトを初期化する。
        
        Args:
            **kwargs: パラメータ
        """
        self.params = kwargs

# インフラ
class REstretto(DockingToolABC):
    """
    REstretto を使ったドッキング計算を行うクラス。
    """
    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        """
        タンパク質について、REstretto用の前処理を行う。
        
        Args:
            protein: 前処理するタンパク質
            
        Returns:
            前処理済みのタンパク質
        """
        raise NotImplementedError()

    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        """
        化合物について、REstretto用の前処理を行う。
        
        Args:
            compound_set: 前処理する化合物セット
            
        Returns:
            前処理済みの化合物セット
        """
        raise NotImplementedError()

    def dock(self, parameters: DockingParameters) -> DockingResult:
        """
        REstrettoを使ってドッキング計算を実施する。
        
        Args:
            parameters: ドッキングパラメータ
            
        Returns:
            ドッキング結果
        """
        raise NotImplementedError()

