"""
docking_automation.molecule パッケージ

分子表現に関するクラスを提供します。
"""

from .protein import Protein
from .compound_set import CompoundSet

__all__ = ["Protein", "CompoundSet"]