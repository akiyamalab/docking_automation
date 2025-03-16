from dataclasses import dataclass, field
from typing import List, Dict, Any, Optional, Tuple

@dataclass(frozen=True)
class Atom:
    """原子を表す値オブジェクト"""
    atom_id: int
    element: str
    x: float
    y: float
    z: float
    charge: float = 0.0
    properties: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class Bond:
    """結合を表す値オブジェクト"""
    atom1_id: int
    atom2_id: int
    bond_type: str  # single, double, triple, aromatic など
    properties: Dict[str, Any] = field(default_factory=dict)


@dataclass(frozen=True)
class MoleculeStructure:
    """分子構造を表す値オブジェクト"""
    atoms: List[Atom]
    bonds: List[Bond]
    properties: Dict[str, Any] = field(default_factory=dict)
    
    def validate(self) -> bool:
        """分子構造の妥当性を検証"""
        # 原子IDが重複していないことを確認
        atom_ids = [atom.atom_id for atom in self.atoms]
        if len(atom_ids) != len(set(atom_ids)):
            return False
        
        # 結合が有効な原子IDを参照していることを確認
        valid_atom_ids = set(atom_ids)
        for bond in self.bonds:
            if bond.atom1_id not in valid_atom_ids or bond.atom2_id not in valid_atom_ids:
                return False
        
        return True
    
    def get_center_of_mass(self) -> Tuple[float, float, float]:
        """分子の重心を計算"""
        if not self.atoms:
            return (0.0, 0.0, 0.0)
        
        total_x = sum(atom.x for atom in self.atoms)
        total_y = sum(atom.y for atom in self.atoms)
        total_z = sum(atom.z for atom in self.atoms)
        n_atoms = len(self.atoms)
        
        return (total_x / n_atoms, total_y / n_atoms, total_z / n_atoms)
    
    def get_bounding_box(self) -> Tuple[float, float, float, float, float, float]:
        """分子の境界ボックスを計算 (min_x, min_y, min_z, max_x, max_y, max_z)"""
        if not self.atoms:
            return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        
        min_x = min(atom.x for atom in self.atoms)
        min_y = min(atom.y for atom in self.atoms)
        min_z = min(atom.z for atom in self.atoms)
        max_x = max(atom.x for atom in self.atoms)
        max_y = max(atom.y for atom in self.atoms)
        max_z = max(atom.z for atom in self.atoms)
        
        return (min_x, min_y, min_z, max_x, max_y, max_z)