from dataclasses import dataclass, field
from datetime import datetime
from typing import ClassVar

from docking_automation.domain.domain_event import DomainEvent

@dataclass(frozen=True)
class ProteinRegistered(DomainEvent):
    """
    タンパク質が登録されたことを表すイベント。
    """
    event_type: ClassVar[str] = "ProteinRegistered"
    protein_id: str = field(default="")
    path: str = field(default="")
    
    def __post_init__(self):
        if not self.protein_id:
            raise ValueError("protein_id は必須です")
        if not self.path:
            raise ValueError("path は必須です")

@dataclass(frozen=True)
class ProteinPathUpdated(DomainEvent):
    """
    タンパク質のパスが更新されたことを表すイベント。
    """
    event_type: ClassVar[str] = "ProteinPathUpdated"
    protein_id: str = field(default="")
    old_path: str = field(default="")
    new_path: str = field(default="")
    
    def __post_init__(self):
        if not self.protein_id:
            raise ValueError("protein_id は必須です")
        if not self.old_path:
            raise ValueError("old_path は必須です")
        if not self.new_path:
            raise ValueError("new_path は必須です")