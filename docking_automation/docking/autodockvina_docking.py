from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from .docking_result import DockingResult
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from .docking_parameters import DockingParameters, SpecificDockingParametersABC
from .docking import DockingToolABC

# 値オブジェクト
class AutoDockVinaParameters(SpecificDockingParametersABC):
    ...

# インフラ
class AutoDockVina(DockingToolABC):
    """
    AutoDock Vina を使ったドッキング計算を行うクラス。
    """
    def preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        ...

    def preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        ...

    def dock(self, parameters: DockingParameters) -> DockingResult:
        ...
