#!/usr/bin/env python
# coding: utf-8

"""
ドッキング結果をSDFファイルに変換するスクリプト

Vinaのドッキング結果（PDBQT形式）をSDFファイルに変換して、
スコア情報を含めた1つのSDFファイルにまとめます。
Open Babelを使用してファイル変換を行います。
"""

import os
import sys
import glob
import tempfile
import subprocess
from pathlib import Path
import re
from rdkit import Chem
from rdkit.Chem import AllChem

# 入出力ディレクトリ
RESULTS_DIR = "output/ALDR/results"
OUTPUT_FILE = "output/ALDR/docking_results.sdf"

def extract_score_from_log(log_file):
    """
    ログファイルからベストスコアを抽出
    
    Args:
        log_file: ログファイルパス
        
    Returns:
        ベストスコア（浮動小数点数）またはNone
    """
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
        
        # "mode |   affinity"行を探す
        mode_line_idx = -1
        for i, line in enumerate(lines):
            if "mode |   affinity" in line:
                mode_line_idx = i
                break
        
        if mode_line_idx >= 0 and mode_line_idx + 3 < len(lines):
            # 最初のスコア行（ベストスコア）を取得
            score_line = lines[mode_line_idx + 3]
            parts = score_line.strip().split()
            if len(parts) >= 2:
                return float(parts[1])
        
        return None
    except Exception as e:
        print(f"警告: スコアの抽出に失敗しました: {e}")
        return None

def parse_ligand_name(filename):
    """
    ファイル名からリガンド名を抽出
    
    Args:
        filename: ファイル名
        
    Returns:
        リガンド名
    """
    # docking_CHEMBL12345_out.pdbqt -> CHEMBL12345
    match = re.search(r'docking_([^_]+)_out\.pdbqt', filename)
    if match:
        return match.group(1)
    return os.path.basename(filename)

def extract_first_model_from_pdbqt(pdbqt_file, output_file=None):
    """
    PDBQTファイルから最初のモデル（最もスコアの良いポーズ）を抽出する
    
    Args:
        pdbqt_file: 入力PDBQTファイルのパス
        output_file: 出力ファイルのパス（省略時は一時ファイル）
        
    Returns:
        抽出したモデルを保存したファイルのパス、変換失敗時はNone
    """
    try:
        # 一時ファイルを作成
        if output_file is None:
            temp_dir = tempfile.mkdtemp()
            output_file = os.path.join(temp_dir, os.path.basename(pdbqt_file))
        
        # PDBQTファイルを読み込む
        with open(pdbqt_file, 'r') as f:
            pdbqt_lines = f.readlines()
        
        # 最初のモデルのみ抽出
        model_lines = []
        in_model = False
        model_complete = False
        
        for line in pdbqt_lines:
            if line.startswith('MODEL') and not in_model and not model_complete:
                in_model = True
                model_lines.append(line)
            elif line.startswith('ENDMDL') and in_model:
                in_model = False
                model_complete = True
                model_lines.append(line)
            elif in_model:
                model_lines.append(line)
            elif not model_complete and not line.startswith('MODEL'):
                # モデルタグがない場合（単一モデル）
                model_lines.append(line)
        
        # ファイルに書き込む
        with open(output_file, 'w') as f:
            f.writelines(model_lines)
        
        return output_file
    except Exception as e:
        print(f"エラー: PDBQTファイルからモデル抽出に失敗しました ({pdbqt_file}): {e}")
        return None

def convert_to_sdf(input_file, output_file=None, input_format='pdbqt'):
    """
    ファイルをSDF形式に変換
    
    Args:
        input_file: 入力ファイルのパス
        output_file: 出力ファイルのパス（省略時は一時ファイル）
        input_format: 入力ファイルの形式（デフォルト: 'pdbqt'）
        
    Returns:
        変換したSDFファイルのパス、変換失敗時はNone
    """
    try:
        if output_file is None:
            temp_dir = tempfile.mkdtemp()
            output_file = os.path.join(temp_dir, os.path.basename(input_file).replace(f'.{input_format}', '.sdf'))
        
        # Open Babelを使用して変換
        cmd = ["obabel", str(input_file), f"-i{input_format}", f"-osdf", "-O", str(output_file)]
        
        print(f"実行コマンド: {' '.join(cmd)}")
        
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        
        if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
            return output_file
        else:
            print(f"エラー: SDFファイルの作成に失敗しました: {output_file}")
            return None
    except Exception as e:
        print(f"エラー: ファイル変換に失敗しました ({input_file}): {e}")
        return None

def main():
    """メイン関数"""
    print(f"=== ドッキング結果をSDFファイルに変換 ===")
    
    # 結果ファイルの検索
    pdbqt_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "docking_*_out.pdbqt")))
    if not pdbqt_files:
        print(f"エラー: {RESULTS_DIR} にPDBQTファイルが見つかりません")
        sys.exit(1)
    
    print(f"検出された結果ファイル: {len(pdbqt_files)}ファイル")
    
    # 結果を格納するリスト
    molecules = []
    scores = {}
    
    # ログファイルからスコアを抽出
    for pdbqt_file in pdbqt_files:
        ligand_name = parse_ligand_name(pdbqt_file)
        log_file = pdbqt_file.replace("_out.pdbqt", "_log.txt")
        
        if os.path.exists(log_file):
            score = extract_score_from_log(log_file)
            if score is not None:
                scores[ligand_name] = score
                print(f"リガンド {ligand_name}: スコア = {score:.2f} kcal/mol")
    
    # スコアでソート
    sorted_scores = sorted(scores.items(), key=lambda x: x[1])
    for i, (name, score) in enumerate(sorted_scores):
        print(f"順位 {i+1}: {name} ({score:.2f} kcal/mol)")
    # 各PDBQTファイルをSDFに変換
    temp_sdf_files = []
    
    for pdbqt_file in pdbqt_files:
        ligand_name = parse_ligand_name(pdbqt_file)
        
        # 最初のモデルを抽出
        print(f"リガンド {ligand_name} を処理中...")
        temp_pdbqt = extract_first_model_from_pdbqt(pdbqt_file)
        
        if temp_pdbqt is not None:
            # SDFに変換
            temp_sdf = convert_to_sdf(temp_pdbqt)
            
            if temp_sdf is not None:
                temp_sdf_files.append((ligand_name, temp_sdf))
                print(f"リガンド {ligand_name} をSDFに変換しました")
            
            # 一時PDBQTファイルを削除
            if os.path.exists(temp_pdbqt):
                os.remove(temp_pdbqt)
            print(f"リガンド {ligand_name} をSDFに変換しました")
    
    # 全SDFファイルを1つにマージして、スコア情報を付加
    writer = Chem.SDWriter(OUTPUT_FILE)
    
    for ligand_name, sdf_file in temp_sdf_files:
        try:
            # SDFファイルを読み込む
            suppl = Chem.SDMolSupplier(sdf_file)
            
            # ベストポーズ（最初のモデル）を取得
            for i, mol in enumerate(suppl):
                if mol is not None and i == 0:  # 最初のモデルのみ処理
                    # プロパティの設定
                    mol.SetProp("_Name", ligand_name)
                    
                    # スコアとランキングをプロパティとして追加
                    if ligand_name in scores:
                        score = scores[ligand_name]
                        mol.SetProp("VINA_SCORE", f"{score:.2f}")
                        
                        # ランキングを追加
                        for i, (name, _) in enumerate(sorted_scores):
                            if name == ligand_name:
                                mol.SetProp("RANK", str(i+1))
                                break
                    
                    # SDFファイルに書き込む
                    writer.write(mol)
                    break
            
            # 一時ファイルを削除
            os.remove(sdf_file)
            
        except Exception as e:
            print(f"警告: SDFファイル処理に失敗しました ({sdf_file}): {e}")
    
    writer.close()
    print(f"\nすべての結果を {OUTPUT_FILE} に保存しました")
    
    # 分子数を確認
    try:
        suppl = Chem.SDMolSupplier(OUTPUT_FILE)
        mol_count = len([mol for mol in suppl if mol is not None])
        print(f"SDFファイルに含まれる分子数: {mol_count}")
    except Exception as e:
        print(f"警告: SDFファイルの確認に失敗しました: {e}")

if __name__ == "__main__":
    main()