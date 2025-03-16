from typing import List, Dict, Tuple, Callable, Any
from sortedcontainers import SortedList

from .docking_result import DockingResult


class DockingResultCollection:
    """
    複数のDockingResultを管理するコレクションクラス。
    
    常にドッキングスコアでソートされた状態を維持する。
    """
    
    def __init__(self):
        """
        DockingResultCollectionオブジェクトを初期化する。
        """
        self._results = SortedList(key=lambda result: result.get_score())
    
    def add(self, result: DockingResult) -> None:
        """
        結果を追加する。
        
        Args:
            result: 追加する結果
        """
        raise NotImplementedError()
    
    def get_all(self) -> List[DockingResult]:
        """
        すべての結果を取得する。
        
        Returns:
            結果のリスト
        """
        raise NotImplementedError()
    
    def get_top(self, n: int) -> List[DockingResult]:
        """
        スコアが上位n件の結果を取得する。
        
        Args:
            n: 取得する件数
            
        Returns:
            上位n件の結果のリスト
        """
        raise NotImplementedError()
    
    def filter(self, condition: Callable[[DockingResult], bool]) -> 'DockingResultCollection':
        """
        条件に合致する結果のみを含む新しいコレクションを作成する。
        
        Args:
            condition: フィルタリング条件
            
        Returns:
            フィルタリングされた結果のコレクション
        """
        raise NotImplementedError()
    
    def group_by_protein(self) -> Dict[str, 'DockingResultCollection']:
        """
        タンパク質IDごとに結果をグループ化する。
        
        Returns:
            タンパク質IDをキー、結果のコレクションを値とする辞書
        """
        raise NotImplementedError()
    
    def group_by_compound(self) -> Dict[Tuple[str, int], 'DockingResultCollection']:
        """
        化合物（化合物セットIDと化合物インデックスの組）ごとに結果をグループ化する。
        
        Returns:
            (化合物セットID, 化合物インデックス)をキー、結果のコレクションを値とする辞書
        """
        raise NotImplementedError()