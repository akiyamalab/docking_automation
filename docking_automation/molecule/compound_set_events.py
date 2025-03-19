from dataclasses import dataclass, field
from datetime import datetime
from typing import ClassVar

from docking_automation.domain.domain_event import DomainEvent

@dataclass(frozen=True)
class CompoundSetRegistered(DomainEvent):
    """
    化合物セットが登録されたことを表すイベント。
    """
    event_type: ClassVar[str] = "CompoundSetRegistered"
    compound_set_id: str = field(default="")
    path: str = field(default="")
    
    def __post_init__(self):
        if not self.compound_set_id:
            raise ValueError("compound_set_id は必須です")
        if not self.path:
            raise ValueError("path は必須です")

@dataclass(frozen=True)
class CompoundSetPathUpdated(DomainEvent):
    """
    化合物セットのパスが更新されたことを表すイベント。
    """
    event_type: ClassVar[str] = "CompoundSetPathUpdated"
    compound_set_id: str = field(default="")
    old_path: str = field(default="")
    new_path: str = field(default="")
    
    def __post_init__(self):
        if not self.compound_set_id:
            raise ValueError("compound_set_id は必須です")
        if not self.old_path:
            raise ValueError("old_path は必須です")
        if not self.new_path:
            raise ValueError("new_path は必須です")