from abc import ABC

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
        protein_file_path, 
        compound_set_file_path,
        grid_box,
    ):
        self.protein = protein
        self.compound_set = compound_set
        self.protein_file_path = protein_file_path
        self.compound_set_file_path = compound_set_file_path
        self.grid_box = grid_box
        
        
    
# 値オブジェクト
class SpecificDockingParametersABC(ABC):
    """
    各ツール固有のパラメータを保持するクラス。
    """
    ...

# 値オブジェクト
class DockingParameters:
    """
    ドッキングツールに関するパラメータを保持するクラス。
    """
    def __init__(
        self, 
        common: CommonDockingParameters, 
        specific: SpecificDockingParametersABC
    ):
        self.common = common
        self.specific = specific