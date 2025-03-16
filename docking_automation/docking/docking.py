
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from .docking_result import DockingResult
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from abc import ABC, abstractmethod

from .docking_parameters import DockingParameters

# インフラ
class DockingToolABC(ABC):
    @abstractmethod
    def preprocess_protein(self, protein: Protein) -> PreprocessedProtein:

        """
        タンパク質について、ドッキング計算の為の前処理を行う。
        """
        ...

    @abstractmethod
    def preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        """
        化合物について、ドッキング計算の為の前処理を行う。
        """
        ...

    @abstractmethod
    def dock(self, parameters: DockingParameters) -> DockingResult:
        """
        ドッキング計算を実施する。
        """
        ...
