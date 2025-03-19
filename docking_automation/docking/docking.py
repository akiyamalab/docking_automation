from typing import List

from docking_automation.docking.preprocessed_compound_set import PreprocessedCompoundSet
from docking_automation.docking.preprocessed_protein import PreprocessedProtein
from .docking_result import DockingResult
from .docking_result_collection import DockingResultCollection
from ..molecule.compound_set import CompoundSet
from ..molecule.protein import Protein
from abc import ABC, abstractmethod
from .grid_box import GridBox

from .docking_parameters import DockingParameters, CommonDockingParameters, SpecificDockingParametersABC

# インフラ
class DockingToolABC(ABC):
    @abstractmethod
    def _preprocess_protein(self, protein: Protein) -> PreprocessedProtein:
        """
        タンパク質について、ドッキング計算の為の前処理を行う。
        
        Args:
            protein: 前処理するタンパク質
            
        Returns:
            前処理済みのタンパク質
        """
        raise NotImplementedError()

    @abstractmethod
    def _preprocess_compound_set(self, compound_set: CompoundSet) -> PreprocessedCompoundSet:
        """
        化合物について、ドッキング計算の為の前処理を行う。
        
        Args:
            compound_set: 前処理する化合物セット
            
        Returns:
            前処理済みの化合物セット
        """
        raise NotImplementedError()

    @abstractmethod
    def dock(self, parameters: DockingParameters) -> List[DockingResult]:
        """
        ドッキング計算を実施する。
        
        複数の化合物に対してドッキング計算を行い、結果のリストを返す。
        
        Args:
            parameters: ドッキングパラメータ
            
        Returns:
            ドッキング結果のリスト
        """
        raise NotImplementedError()
    
    def _dock_with_auto_preprocess(
        self,
        protein: Protein,
        compound_set: CompoundSet,
        grid_box: GridBox,
        additional_params: SpecificDockingParametersABC
    ) -> List[DockingResult]:
        """
        前処理を自動的に行ってからドッキング計算を実施する。
        内部メソッドであり、直接呼び出すべきではない。
        
        複数の化合物に対してドッキング計算を行い、結果のリストを返す。
        
        Args:
            protein: ドッキング対象のタンパク質
            compound_set: ドッキング対象の化合物セット
            grid_box: ドッキング計算の領域
            additional_params: 追加のパラメータ
            
        Returns:
            ドッキング結果のリスト
        """
        # 前処理
        preprocessed_protein = self._preprocess_protein(protein)
        preprocessed_compound_set = self._preprocess_compound_set(compound_set)
        # パラメータの設定
        common_params = CommonDockingParameters(
            protein=preprocessed_protein,
            compound_set=preprocessed_compound_set,
            grid_box=grid_box
        )
        
        # 具体的なパラメータの設定は各ドッキングツールの実装に依存
        parameters = DockingParameters(
            common=common_params,
            specific=additional_params
        )
        
        # ドッキング実行
        return self.dock(parameters)
    
    def run_docking(
        self,
        protein: Protein,
        compound_set: CompoundSet,
        grid_box: GridBox,
        additional_params: SpecificDockingParametersABC
    ) -> DockingResultCollection:
        """
        ドッキング計算を実行し、結果をコレクションとして返す。
        
        複数の化合物に対してドッキング計算を行い、結果をコレクションに追加する。
        
        Args:
            protein: ドッキング対象のタンパク質
            compound_set: ドッキング対象の化合物セット
            grid_box: ドッキング計算の領域
            additional_params: 追加のパラメータ
            
        Returns:
            ドッキング結果のコレクション
        """
        results = self._dock_with_auto_preprocess(
            protein=protein,
            compound_set=compound_set,
            grid_box=grid_box,
            additional_params=additional_params
        )
        
        # 結果をコレクションに変換
        collection = DockingResultCollection()
        collection.extend(results)  # extendを使用して複数の結果を一度に追加
        
        return collection
