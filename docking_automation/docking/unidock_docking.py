from __future__ import annotations

from typing import List

from docking_automation.docking.docking import DockingToolABC
from docking_automation.docking.docking_parameters import DockingParameters, UniDockParameters
from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection
from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein
from docking_automation.docking.grid_box import GridBox


class UniDockDocking(DockingToolABC):
    """Uni-Dock GPU batch docking backend."""

    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        raise NotImplementedError("UniDockDocking._preprocess_protein: implemented in 011_b")

    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        raise NotImplementedError("UniDockDocking._preprocess_compound_set: implemented in 011_b")

    def dock(self, parameters: DockingParameters) -> List[DockingResult]:
        raise NotImplementedError("UniDockDocking.dock: implemented in 011_b")

    def run_docking(
        self,
        protein: Protein,
        compound_set: CompoundSet,
        grid_box: GridBox,
        parameters: UniDockParameters | None = None,
    ) -> DockingResultCollection:
        """バッチ GPU ドッキング実行。"""
        if parameters is None:
            parameters = UniDockParameters()
        raise NotImplementedError("UniDockDocking.run_docking: implemented in 011_b")
