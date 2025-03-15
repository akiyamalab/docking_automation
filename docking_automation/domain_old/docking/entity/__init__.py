"""ドッキングドメインのエンティティモジュール

このモジュールはドッキング計算に関するエンティティを提供します。
"""

# ドッキング関連のエンティティをエクスポート
from docking_automation.domain.docking.entity.docking_result import DockingResult
from docking_automation.domain.docking.entity.docking_task import DockingTask
from docking_automation.domain.docking.entity.ligand import Ligand
from docking_automation.domain.docking.entity.receptor import Receptor

__all__ = [
    'DockingResult',
    'DockingTask',
    'Ligand',
    'Receptor',
]