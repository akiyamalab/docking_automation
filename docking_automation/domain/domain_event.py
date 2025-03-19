from dataclasses import dataclass, field
from datetime import datetime
from typing import ClassVar, Dict, Any

@dataclass(frozen=True)
class DomainEvent:
    """
    ドメインイベントの基底クラス。
    
    ドメインで発生した重要な出来事を表現するイベントオブジェクト。
    不変（イミュータブル）であり、過去に発生した事実を表現する。
    """
    occurred_on: datetime = field(default_factory=datetime.now)
    event_type: ClassVar[str] = "DomainEvent"
    
    def to_dict(self) -> Dict[str, Any]:
        """
        イベントを辞書形式に変換する。
        
        Returns:
            イベントの内容を表す辞書
        """
        result = {
            "event_type": self.event_type,
            "occurred_on": self.occurred_on.isoformat()
        }
        
        # dataclassの全フィールドを辞書に追加
        for field_name, field_value in self.__dict__.items():
            if field_name != "occurred_on":  # 既に追加済みのフィールドはスキップ
                result[field_name] = field_value
                
        return result