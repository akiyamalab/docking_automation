"""ドッキングドメインのエンティティモジュール

このモジュールはドッキング計算に関するエンティティを提供します。
"""

from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.entity.docking_task import DockingTask, TaskStatus
from docking_automation.docking.entity.docking_result import DockingResult

__all__ = [
    'Ligand',
    'Receptor',
    'DockingTask',
    'TaskStatus',
    'DockingResult',
]