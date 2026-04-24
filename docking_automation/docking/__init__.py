"""
docking_automation.docking パッケージ

ドッキング計算に関するクラスを提供します。

各ツール (AutoDockVina / UniDockDocking / UniDock2Docking) は lazy import で
公開する。これにより、例えば UniDock2Docking を使う環境で openbabel が
無くてもインポートエラーにならない。また、openbabel と msys の C++ stdlib
ABI 衝突で segfault する問題 (unidock2 env で観察) を回避できる。
"""
from __future__ import annotations

from typing import TYPE_CHECKING

# 軽量な値オブジェクト群は eager import で OK (外部 C ライブラリ依存なし)
from .docking_parameters import CommonDockingParameters, DockingParameters
from .docking_result import DockingResult
from .docking_result_collection import DockingResultCollection
from .grid_box import GridBox

if TYPE_CHECKING:
    from .autodockvina_docking import AutoDockVina, AutoDockVinaParameters


def __getattr__(name: str):
    """AutoDockVina など外部依存のあるシンボルは lazy import。"""
    if name in ('AutoDockVina', 'AutoDockVinaParameters'):
        from .autodockvina_docking import AutoDockVina, AutoDockVinaParameters
        return {
            'AutoDockVina': AutoDockVina,
            'AutoDockVinaParameters': AutoDockVinaParameters,
        }[name]
    raise AttributeError(f'module {__name__!r} has no attribute {name!r}')


__all__ = [
    'GridBox',
    'AutoDockVina',
    'AutoDockVinaParameters',
    'DockingResult',
    'DockingResultCollection',
    'DockingParameters',
    'CommonDockingParameters',
]
