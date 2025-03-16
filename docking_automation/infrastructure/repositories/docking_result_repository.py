from typing import List, Tuple

from docking_automation.docking.docking_result import DockingResult
from docking_automation.docking.docking_result_collection import DockingResultCollection


# リポジトリ
class DockingResultRepository:
    """ドッキング計算結果を保持するリポジトリクラス。"""
    
    def save(self, result: DockingResult) -> None:
        """
        ドッキング計算結果を保存する。
        
        Args:
            result: 保存する結果
        """
        raise NotImplementedError()
    
    def load(self, result_id: str) -> DockingResult:
        """
        ドッキング計算結果を読み込む。
        
        Args:
            result_id: 結果のID
            
        Returns:
            読み込んだ結果
        """
        raise NotImplementedError()
    
    def find_by_protein(self, protein_id: str) -> DockingResultCollection:
        """
        タンパク質IDに基づいて結果を検索する。
        
        Args:
            protein_id: タンパク質のID
            
        Returns:
            検索結果のコレクション
        """
        raise NotImplementedError()
    
    def find_by_compound_set(self, compound_set_id: str) -> DockingResultCollection:
        """
        化合物セットIDに基づいて結果を検索する。
        
        Args:
            compound_set_id: 化合物セットのID
            
        Returns:
            検索結果のコレクション
        """
        raise NotImplementedError()
    
    def find_by_compound(self, compound_set_id: str, compound_index: int) -> DockingResultCollection:
        """
        化合物（化合物セットIDと化合物インデックスの組）に基づいて結果を検索する。
        
        Args:
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            
        Returns:
            検索結果のコレクション
        """
        raise NotImplementedError()
    
    def find_by_protein_and_compound(self, protein_id: str, compound_set_id: str, compound_index: int) -> DockingResultCollection:
        """
        タンパク質IDと化合物（化合物セットIDと化合物インデックスの組）に基づいて結果を検索する。
        
        Args:
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            
        Returns:
            検索結果のコレクション
        """
        raise NotImplementedError()
