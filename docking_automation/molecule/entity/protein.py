"""タンパク質エンティティ

このモジュールは、タンパク質を表すエンティティクラスを提供します。
"""

import os
from typing import Dict, Any, Optional, Set


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
    
    def __str__(self) -> str:
        """文字列表現を返す
        
        Returns:
            オブジェクトの文字列表現
        """
        return f"Protein({self.id}, name={self.name}, chains={len(self.chains)})"