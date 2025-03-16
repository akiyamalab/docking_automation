"""タンパク質エンティティ

このモジュールは、タンパク質を表すエンティティクラスを提供します。
"""

import os
from typing import Dict, Any, Optional, Set, List, Union, Tuple

from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom
from docking_automation.docking.value_object.grid_box import GridBox


class Protein:
    """タンパク質エンティティ
    
    ドッキング計算で使用されるタンパク質（レセプター）を表すエンティティです。
    タンパク質の構造情報やアクティブサイト情報を保持します。
    """
    
    def __init__(
        self,
        id: str,
        path: str,
        structure: Optional[Any] = None,
        format: Optional[Any] = None,
        chains: Optional[Set[str]] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """コンストラクタ
        
        Args:
            id: タンパク質の一意識別子
            path: タンパク質ファイルのパス
            structure: タンパク質の構造情報（オプション）
            format: タンパク質のファイル形式情報（オプション）
            chains: タンパク質のチェーン識別子（オプション）
            metadata: 追加情報（辞書形式）
        """
        self.id = id
        self.path = path
        self.structure = structure
        self.format = format
        self.chains = chains or set()
        self.metadata = metadata or {}
        # ファイル名をデフォルト名として使用
        self._name = os.path.splitext(os.path.basename(path))[0]
    
    @property
    def name(self) -> str:
        """タンパク質の名前を取得する
        
        Returns:
            タンパク質の名前
        """
        return self._name
    
    @name.setter
    def name(self, value: str) -> None:
        """タンパク質の名前を設定する
        
        Args:
            value: 設定する名前
        """
        self._name = value
    
    def add_metadata(self, metadata: Dict[str, Any]) -> None:
        """メタデータを追加する
        
        Args:
            metadata: 追加するメタデータ辞書
        """
        self.metadata.update(metadata)
    
    def has_active_site_defined(self) -> bool:
        """アクティブサイトが定義されているかどうかを確認する
        
        Returns:
            アクティブサイトが定義されている場合はTrue、そうでない場合はFalse
        """
        # アクティブサイトが定義されているかどうかをメタデータなどから判定
        return 'active_site' in self.metadata
    
    def get_active_site_residues(self) -> List[str]:
        """アクティブサイトの残基IDリストを取得する
        
        Returns:
            アクティブサイト残基IDのリスト（定義されていない場合は空リスト）
        """
        if not self.has_active_site_defined():
            return []
        
        # メタデータからアクティブサイト情報を取得
        active_site = self.metadata.get('active_site', {})
        return active_site.get('residues', [])
    
    def get_residue_atoms(self, residue_ids: List[str]) -> Optional[MoleculeStructure]:
        """指定された残基IDに対応する原子を含む構造を取得する
        
        Args:
            residue_ids: 対象となる残基IDのリスト
        
        Returns:
            指定された残基の原子を含む構造オブジェクト（見つからない場合はNone）
        
        Notes:
            現状では簡易実装のため、構造データがない場合はNoneを返します。
            将来的には、構造データから残基IDを解析して原子を抽出する処理が必要です。
        """
        # 構造データがない場合はNoneを返す
        if not hasattr(self, 'structure') or self.structure is None:
            return None
        
        # TODO: PDB/CIFパーサーを実装し、残基IDから対応する原子を抽出する処理を追加
        # 現状では単純に元の構造をそのまま返す（仮実装）
        return self.structure
    
    def calculate_grid_box(
        self,
        active_site_only: bool = False,
        padding: float = 4.0,
        default_box_size: Optional[float] = 20.0
    ) -> GridBox:
        """タンパク質からグリッドボックスを計算する
        
        Args:
            active_site_only: アクティブサイト情報のみを使用するかどうか
            padding: 境界ボックスに追加するパディング（Å）
            default_box_size: デフォルトのボックスサイズ（Å）
        
        Returns:
            計算されたグリッドボックス
        
        Raises:
            ValueError: 構造データがない場合
            
        Examples:
            >>> protein = Protein(id="1abc", path="/path/to/1abc.pdb")
            >>> grid_box = protein.calculate_grid_box()
            >>> print(f"中心座標: {grid_box.get_center()}")
            >>> print(f"ボックスサイズ: {grid_box.get_size()}")
        """
        from docking_automation.molecule.service.grid_box_calculator import GridBoxCalculator
        
        # アクティブサイトのみを使用する場合
        if active_site_only and self.has_active_site_defined():
            residue_ids = self.get_active_site_residues()
            return GridBoxCalculator.calculate_from_protein(
                self,
                active_site_residues=residue_ids,
                padding=padding,
                default_box_size=None  # アクティブサイトの場合は実際のサイズを使用
            )
        
        # 通常のケース：タンパク質全体から計算
        return GridBoxCalculator.calculate_from_protein(
            self,
            padding=padding,
            default_box_size=default_box_size
        )
    
    def __str__(self) -> str:
        """文字列表現を返す
        
        Returns:
            オブジェクトの文字列表現
        """
        return f"Protein({self.id}, name={self.name}, chains={len(self.chains)})"