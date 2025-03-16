#!/usr/bin/env python
# coding: utf-8

"""
分子データからドッキング計算のグリッドボックスを計算するサンプルスクリプト

様々な分子形式（PDB、MOL2など）からグリッドボックスを計算し、結果を表示します。
"""

import os
import sys
import argparse
import numpy as np
from typing import Optional, Tuple, List, Dict, Any

# 必要なモジュールのインポート
from docking_automation.molecule.entity.protein import Protein
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.molecule.service.grid_box_calculator import GridBoxCalculator
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond


def parse_mol2(filename: str) -> MoleculeStructure:
    """MOL2ファイルから分子構造を作成する
    
    Args:
        filename: MOL2ファイルのパス
        
    Returns:
        MoleculeStructure: 分子構造オブジェクト
    """
    atoms = []
    bonds = []
    atom_section = False
    bond_section = False
    atom_id_map = {}  # 元のID → 新しいID
    
    with open(filename, 'r') as f:
        lines = f.readlines()
        
        for i, line in enumerate(lines):
            line = line.strip()
            
            # セクションの判定
            if line == '@<TRIPOS>ATOM':
                atom_section = True
                bond_section = False
                continue
            elif line == '@<TRIPOS>BOND':
                atom_section = False
                bond_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                atom_section = False
                bond_section = False
                continue
            
            # 原子セクションの処理
            if atom_section and line:
                parts = line.split()
                if len(parts) >= 9:
                    atom_id = int(parts[0])
                    atom_name = parts[1]
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atom_type = parts[5]
                    
                    # Charge情報の取得（あれば）
                    charge = 0.0
                    if len(parts) > 8:
                        try:
                            charge = float(parts[8])
                        except ValueError:
                            pass
                    
                    # 新しいID（1-based）
                    new_id = len(atoms) + 1
                    atom_id_map[atom_id] = new_id
                    
                    # 原子オブジェクトの作成
                    atoms.append(Atom(
                        atom_id=new_id,
                        element=atom_type.split('.')[0] if '.' in atom_type else atom_type,
                        x=x, y=y, z=z,
                        charge=charge,
                        properties={'name': atom_name, 'type': atom_type}
                    ))
            
            # 結合セクションの処理
            if bond_section and line:
                parts = line.split()
                if len(parts) >= 4:
                    bond_id = int(parts[0])
                    atom1_id = int(parts[1])
                    atom2_id = int(parts[2])
                    bond_type = parts[3]
                    
                    # ID変換
                    atom1_id_new = atom_id_map.get(atom1_id, atom1_id)
                    atom2_id_new = atom_id_map.get(atom2_id, atom2_id)
                    
                    # 結合タイプの変換
                    bond_type_map = {
                        '1': 'single',
                        '2': 'double',
                        '3': 'triple',
                        'ar': 'aromatic',
                        'am': 'amide',
                        'du': 'dummy',
                        'un': 'unknown',
                        'nc': 'not connected'
                    }
                    bond_type_str = bond_type_map.get(bond_type.lower(), 'single')
                    
                    # 結合オブジェクトの作成
                    bonds.append(Bond(
                        atom1_id=atom1_id_new,
                        atom2_id=atom2_id_new,
                        bond_type=bond_type_str
                    ))
    
    # 分子構造の作成
    if not atoms:
        raise ValueError(f"エラー: {filename}から原子を読み取れませんでした")
    
    return MoleculeStructure(atoms=atoms, bonds=bonds)


def create_mock_protein(structure: MoleculeStructure, protein_id: str = "mock") -> Protein:
    """分子構造からモック用のタンパク質オブジェクトを作成する
    
    Args:
        structure: 分子構造
        protein_id: タンパク質ID
        
    Returns:
        Protein: タンパク質オブジェクト
    """
    protein = Protein(
        id=protein_id,
        path="mock_path.pdb",  # ダミーパス
        chains=set()
    )
    
    # 構造を設定
    protein.structure = structure
    
    return protein


def print_grid_box_info(grid_box: GridBox) -> None:
    """グリッドボックスの情報を表示する
    
    Args:
        grid_box: グリッドボックスオブジェクト
    """
    center = grid_box.get_center()
    size = grid_box.get_size()
    min_corner = grid_box.get_corner_min()
    max_corner = grid_box.get_corner_max()
    volume = grid_box.get_volume()
    
    print(f"グリッドボックス情報:")
    print(f"  中心座標:   ({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f})")
    print(f"  ボックスサイズ: {size[0]:.3f} x {size[1]:.3f} x {size[2]:.3f} Å")
    print(f"  最小角:     ({min_corner[0]:.3f}, {min_corner[1]:.3f}, {min_corner[2]:.3f})")
    print(f"  最大角:     ({max_corner[0]:.3f}, {max_corner[1]:.3f}, {max_corner[2]:.3f})")
    print(f"  体積:       {volume:.3f} Å³")


def calculate_from_mol2(mol2_file: str, padding: float = 4.0, fixed_size: Optional[float] = None) -> GridBox:
    """MOL2ファイルからグリッドボックスを計算する
    
    Args:
        mol2_file: MOL2ファイルパス
        padding: ボックスのパディング（Å）
        fixed_size: 固定サイズ（指定すれば全方向同じサイズに設定）
        
    Returns:
        GridBox: 計算されたグリッドボックス
    """
    print(f"MOL2ファイル: {mol2_file} からグリッドボックスを計算します")
    
    # MOL2ファイルを解析
    structure = parse_mol2(mol2_file)
    print(f"読み込んだ原子数: {len(structure.atoms)}")
    
    # グリッドボックスを計算
    grid_box = GridBoxCalculator.calculate_from_structure(
        structure, 
        padding=padding,
        default_box_size=fixed_size
    )
    
    return grid_box


def main() -> None:
    """メイン関数"""
    parser = argparse.ArgumentParser(description='分子データからドッキング用グリッドボックスを計算')
    parser.add_argument('input_file', help='入力分子ファイル（MOL2形式）')
    parser.add_argument('-p', '--padding', type=float, default=4.0, help='境界ボックスに追加するパディング（Å）')
    parser.add_argument('-s', '--size', type=float, help='固定サイズを使用（全方向同じサイズに）')
    parser.add_argument('-o', '--output', help='出力ファイル（指定すると結果をJSONで保存）')
    
    args = parser.parse_args()
    
    # ファイル存在チェック
    if not os.path.exists(args.input_file):
        print(f"エラー: ファイル {args.input_file} が見つかりません")
        sys.exit(1)
    
    try:
        # グリッドボックスを計算
        grid_box = calculate_from_mol2(
            args.input_file, 
            padding=args.padding,
            fixed_size=args.size
        )
        
        # 結果を表示
        print_grid_box_info(grid_box)
        
        # 出力ファイルが指定されていれば保存
        if args.output:
            import json
            with open(args.output, 'w') as f:
                json.dump(grid_box.to_dict(), f, indent=2)
                print(f"グリッドボックス情報を {args.output} に保存しました")
        
    except Exception as e:
        print(f"エラー: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()