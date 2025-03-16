#!/usr/bin/env python
# coding: utf-8

"""
crystal_ligand.mol2からドッキング中心の座標を計算する
"""

import sys
import numpy as np

def parse_mol2(filename):
    """MOL2ファイルから原子座標を抽出する"""
    atoms = []
    atom_section = False
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            
            if line == '@<TRIPOS>ATOM':
                atom_section = True
                continue
            elif line.startswith('@<TRIPOS>'):
                atom_section = False
                continue
            
            if atom_section and line:
                parts = line.split()
                if len(parts) >= 9:
                    # インデックス, 原子名, x, y, z, 原子タイプ, ...
                    x, y, z = float(parts[2]), float(parts[3]), float(parts[4])
                    atoms.append((x, y, z))
    
    return np.array(atoms)

def calculate_center(coords):
    """原子座標の中心を計算する"""
    return np.mean(coords, axis=0)

def main():
    if len(sys.argv) != 2:
        print(f"使用法: {sys.argv[0]} <mol2ファイル>")
        sys.exit(1)
    
    mol2_file = sys.argv[1]
    coords = parse_mol2(mol2_file)
    
    if len(coords) == 0:
        print(f"エラー: {mol2_file}から原子座標を読み取れませんでした")
        sys.exit(1)
    
    center = calculate_center(coords)
    print(f"結晶リガンドの原子数: {len(coords)}")
    print(f"ドッキング中心の座標: {center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}")
    
    # ボックスサイズを標準で20Åとする
    box_size = 20.0
    print(f"推奨ボックスサイズ: {box_size} x {box_size} x {box_size}")
    
    return center, box_size

if __name__ == "__main__":
    main()