"""化合物エンティティ

このモジュールは、化合物を表すエンティティクラスを提供します。
"""

import os
from typing import Dict, Any, Optional


class Compound:
    """化合物エンティティ
    
    ドッキング計算で使用される化合物（リガンド）を表すエンティティです。
    化合物の構造情報や物理化学的特性を保持します。
    """
    
    def __init__(
        self,
        id: str,
        path: str,
        structure: Optional[Any] = None,
        format: Optional[Any] = None,
        metadata: Optional[Dict[str, Any]] = None
    ):
        """コンストラクタ
        
        Args:
            id: 化合物の一意識別子
            path: 化合物ファイルのパス
            structure: 化合物の構造情報（オプション）
            format: 化合物のファイル形式情報（オプション）
            metadata: 追加情報（辞書形式）
        """
        self.id = id
        self.path = path
        self.structure = structure
        self.format = format
        self.metadata = metadata or {}
        # ファイル名をデフォルト名として使用
        self._name = os.path.splitext(os.path.basename(path))[0]
    
    @property
    def name(self) -> str:
        """化合物の名前を取得する
        
        Returns:
            化合物の名前
        """
        return self._name
    
    @name.setter
    def name(self, value: str) -> None:
        """化合物の名前を設定する
        
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
    
    def __str__(self) -> str:
        """文字列表現を返す
        
        Returns:
            オブジェクトの文字列表現
        """
        return f"Compound({self.id}, name={self.name})"