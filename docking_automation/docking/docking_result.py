from pathlib import Path
from typing import Dict, Any, Optional


# エンティティ
class DockingResult:
    """
    ドッキング計算の結果を保持するクラス。
    1つのタンパク質と1つの化合物のドッキング計算結果が保持される。
    """
    
    def __init__(
        self,
        result_path: Path,
        protein_id: str,
        compound_set_id: str,
        compound_index: int,
        docking_score: float,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """
        DockingResultオブジェクトを初期化する。
        
        Args:
            result_path: 結果 SDF ファイルのパス
            protein_id: タンパク質のID
            compound_set_id: 化合物セットのID
            compound_index: 化合物セット内の化合物のインデックス
            docking_score: ドッキングスコア
            metadata: メタデータ
        """
        self.result_path = result_path
        self.protein_id = protein_id
        self.compound_set_id = compound_set_id
        self.compound_index = compound_index
        self.docking_score = docking_score
        self.metadata = metadata or {}