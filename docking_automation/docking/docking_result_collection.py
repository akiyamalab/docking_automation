from typing import List, Dict, Tuple, Callable, Iterator
from sortedcontainers import SortedList # type: ignore
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
        self._results = SortedList(key=lambda result: result.docking_score)
    
    def add(self, result: DockingResult) -> None:
        """
        結果を追加する。
        
        Args:
            result: 追加する結果
        """
        self._results.add(result)
    
    def get_all(self) -> List[DockingResult]:
        """
        すべての結果を取得する。
        
        Returns:
            結果のリスト
        """
        return list(self._results)
    
    def get_top(self, n: int) -> List[DockingResult]:
        """
        スコアが上位n件の結果を取得する。
        
        Args:
            n: 取得する件数
            
        Returns:
            上位n件の結果のリスト
        """
        # 明示的にリストを作成して型を保証する
        top_n = []
        for i, result in enumerate(self._results):
            if i >= n:
                break
            top_n.append(result)
        return top_n
    
    def filter(self, condition: Callable[[DockingResult], bool]) -> 'DockingResultCollection':
        """
        条件に合致する結果のみを含む新しいコレクションを作成する。
        
        Args:
            condition: フィルタリング条件
            
        Returns:
            フィルタリングされた結果のコレクション
        """
        filtered_collection = DockingResultCollection()
        for result in self._results:
            if condition(result):
                filtered_collection.add(result)
        return filtered_collection
    
    def group_by_protein(self) -> Dict[str, 'DockingResultCollection']:
        """
        タンパク質IDごとに結果をグループ化する。
        
        Returns:
            タンパク質IDをキー、結果のコレクションを値とする辞書
        """
        grouped: Dict[str, DockingResultCollection] = {}
        for result in self._results:
            if result.protein_id not in grouped:
                grouped[result.protein_id] = DockingResultCollection()
            grouped[result.protein_id].add(result)
        return grouped
    
    def group_by_compound(self) -> Dict[Tuple[str, int], 'DockingResultCollection']:
        """
        化合物（化合物セットIDと化合物インデックスの組）ごとに結果をグループ化する。
        
        Returns:
            (化合物セットID, 化合物インデックス)をキー、結果のコレクションを値とする辞書
        """
        grouped: Dict[Tuple[str, int], DockingResultCollection] = {}
        for result in self._results:
            key = (result.compound_set_id, result.compound_index)
            if key not in grouped:
                grouped[key] = DockingResultCollection()
            grouped[key].add(result)
        return grouped
    
    def __len__(self) -> int:
        """
        結果の数を取得する。
        
        Returns:
            結果の数
        """
        return len(self._results)
    
    def __iter__(self) -> Iterator[DockingResult]:
        """
        結果のイテレータを取得する。
        
        Returns:
            結果のイテレータ
        """
        return iter(self._results)