#!/usr/bin/env python
# coding: utf-8

"""
PDBQTファイルの形式を修正するスクリプト

このスクリプトは、Open Babelで生成されたPDBQTファイルを
Vinaが認識できる正しい形式に修正します。
"""

import os
import sys
import glob
from pathlib import Path

def fix_pdbqt_file(input_file, output_file=None):
    """
    PDBQTファイルの形式を修正
    
    Args:
        input_file: 入力PDBQTファイルパス
        output_file: 出力PDBQTファイルパス（省略時は入力ファイルを上書き）
    
    Returns:
        修正したファイルパス
    """
    if output_file is None:
        output_file = input_file
    
    try:
        with open(input_file, 'r') as f:
            lines = f.readlines()
        
        # 既存の行からROOT、BRANCH、ENDROOT、TORSDOFなどの行を探す
        has_root = any(line.strip() == 'ROOT' for line in lines)
        has_endroot = any(line.strip() == 'ENDROOT' for line in lines)
        has_torsdof = any(line.startswith('TORSDOF') for line in lines)
        
        # 原子行とブランチ行を分類
        atom_lines = []
        remark_lines = []
        branch_sections = []
        
        current_section = []
        in_branch = False
        
        for line in lines:
            if line.startswith('REMARK'):
                remark_lines.append(line)
            elif line.startswith('ROOT'):
                pass  # ROOTタグは後で追加するので無視
            elif line.startswith('ENDROOT'):
                pass  # ENDROOTタグは後で追加するので無視
            elif line.startswith('BRANCH'):
                if current_section:
                    atom_lines.extend(current_section)
                    current_section = []
                in_branch = True
                current_section.append(line)
            elif line.startswith('ENDBRANCH'):
                current_section.append(line)
                branch_sections.append(current_section)
                current_section = []
                in_branch = False
            elif line.startswith('TORSDOF'):
                pass  # TORSDOFタグは後で追加するので無視
            elif line.startswith('ATOM') or line.startswith('HETATM'):
                if not in_branch:
                    atom_lines.append(line)
                else:
                    current_section.append(line)
        
        # 残りの部分があれば追加
        if current_section:
            if in_branch:
                branch_sections.append(current_section)
            else:
                atom_lines.extend(current_section)
        
        # 回転可能な結合の数を数える
        branch_count = len(branch_sections)
        
        # 修正したファイルを書き出す
        with open(output_file, 'w') as f:
            # REMARKを先頭に出力
            for line in remark_lines:
                f.write(line)
            
            # ROOT
            f.write('ROOT\n')
            
            # 原子行
            for line in atom_lines:
                f.write(line)
            
            # ENDROOT
            f.write('ENDROOT\n')
            
            # ブランチセクション
            for section in branch_sections:
                for line in section:
                    f.write(line)
            
            # TORSDOF
            f.write(f'TORSDOF {max(1, branch_count)}\n')
        
        print(f"PDBQTファイルを修正しました: {output_file}")
        return output_file
    
    except Exception as e:
        print(f"エラー: PDBQTファイルの修正に失敗しました: {e}")
        return None

def main():
    """メイン関数"""
    if len(sys.argv) < 2:
        print("使用法: python fix_pdbqt_format.py <input_pattern>")
        print("例: python fix_pdbqt_format.py output/ALDR/ligand_*.pdbqt")
        sys.exit(1)
    
    pattern = sys.argv[1]
    files = glob.glob(pattern)
    
    if not files:
        print(f"エラー: パターン {pattern} に一致するファイルが見つかりません")
        sys.exit(1)
    
    print(f"修正対象のファイル数: {len(files)}")
    
    success_count = 0
    for file_path in files:
        result = fix_pdbqt_file(file_path)
        if result:
            success_count += 1
    
    print(f"\n処理完了: {success_count}/{len(files)} ファイルを修正しました")

if __name__ == "__main__":
    main()