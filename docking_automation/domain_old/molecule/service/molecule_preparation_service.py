from abc import ABC, abstractmethod
from typing import List, Optional, Dict, Any

from docking_automation.domain.molecule.entity.compound import Compound
from docking_automation.domain.molecule.entity.protein import Protein
from docking_automation.domain.molecule.value_object.molecule_property import MoleculeProperty


class MoleculePreparationService(ABC):
    """分子準備サービスのインターフェース
    
    このサービスは、ドッキング計算のために化合物とタンパク質を準備することを担当します。
    """
    
    @abstractmethod
    def prepare_ligand(self, compound: Compound, method: Optional[str] = None) -> Compound:
        """化合物をリガンドとして準備する
        
        典型的な準備手順：
        1. 3D構造の生成／最適化
        2. 水素原子の追加
        3. 電荷の計算と割り当て
        4. PDBQTフォーマットへの変換
        
        Args:
            compound: 準備する化合物
            method: 準備方法（"gasteiger", "am1bcc"など、実装によって異なる）
            
        Returns:
            準備された化合物
            
        Raises:
            ValueError: 準備中に無効なパラメータが検出された場合
            RuntimeError: 準備処理に失敗した場合
        """
        pass
    
    @abstractmethod
    def prepare_receptor(self, protein: Protein, method: Optional[str] = None) -> Protein:
        """タンパク質をレセプターとして準備する
        
        典型的な準備手順：
        1. 水分子の削除（オプション）
        2. 水素原子の追加
        3. 電荷の計算と割り当て
        4. PDBQTフォーマットへの変換
        
        Args:
            protein: 準備するタンパク質
            method: 準備方法（実装によって異なる）
            
        Returns:
            準備されたタンパク質
            
        Raises:
            ValueError: 準備中に無効なパラメータが検出された場合
            RuntimeError: 準備処理に失敗した場合
        """
        pass
    
    @abstractmethod
    def calculate_properties(self, compound: Compound) -> List[MoleculeProperty]:
        """化合物の物理化学的特性を計算する
        
        典型的な特性：
        - 分子量
        - LogP
        - 水素結合ドナー/アクセプター数
        - 回転可能な結合の数
        - 極性表面積
        
        Args:
            compound: 特性を計算する化合物
            
        Returns:
            計算された特性のリスト
            
        Raises:
            ValueError: 計算に失敗した場合
        """
        pass
    
    @abstractmethod
    def add_hydrogens(self, molecule: Compound | Protein, ph: float = 7.4) -> Compound | Protein:
        """分子に水素原子を追加する
        
        Args:
            molecule: 水素原子を追加する分子
            ph: 溶液のpH値（デフォルト: 7.4）
            
        Returns:
            水素原子が追加された分子
            
        Raises:
            ValueError: 無効なパラメータが指定された場合
            RuntimeError: 処理に失敗した場合
        """
        pass
    
    @abstractmethod
    def remove_water(self, protein: Protein) -> Protein:
        """タンパク質から水分子を除去する
        
        Args:
            protein: 水分子を除去するタンパク質
            
        Returns:
            水分子が除去されたタンパク質
            
        Raises:
            RuntimeError: 処理に失敗した場合
        """
        pass
    
    @abstractmethod
    def get_supported_methods(self) -> Dict[str, str]:
        """このサービスがサポートする準備方法のリストを取得する
        
        Returns:
            サポートされている準備方法の辞書 {method_name: description}
        """
        pass
