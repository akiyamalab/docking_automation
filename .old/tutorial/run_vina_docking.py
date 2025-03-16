#!/usr/bin/env python
# coding: utf-8

"""
AutoDock Vinaを使用してドッキング計算を実行するスクリプト

このスクリプトは、変換済みのPDBQTファイルを使用して
AutoDock Vinaによるドッキング計算を行います。
"""

import os
import sys
import time
import subprocess
import glob
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt

# 設定
RECEPTOR_PDBQT = "output/ALDR/receptor.pdbqt"
LIGANDS_DIR = "output/ALDR"
RESULTS_DIR = "output/ALDR/results"
VINA_EXECUTABLE = "/advina/1.1.2/bin/vina"
MAX_LIGANDS = 5  # 処理する最大リガンド数

# ドッキング設定
DOCKING_CENTER = (16.611, -7.032, 14.448)  # 結晶リガンドから計算した中心座標
BOX_SIZE = (20.0, 20.0, 20.0)  # ボックスサイズ (x, y, z)
EXHAUSTIVENESS = 8  # 探索の徹底度
NUM_MODES = 9      # 出力するポーズ数
CPU = 1            # 使用するCPU数

# 出力ディレクトリの作成
os.makedirs(RESULTS_DIR, exist_ok=True)

def run_vina_docking(receptor_file, ligand_file, output_prefix):
    """
    AutoDock Vinaを使用してドッキング計算を実行
    
    Args:
        receptor_file: レセプター（タンパク質）のPDBQTファイルパス
        ligand_file: リガンド（化合物）のPDBQTファイルパス
        output_prefix: 出力ファイルの接頭辞
        
    Returns:
        結果ファイルのパスとベストスコアのタプル
    """
    ligand_name = Path(ligand_file).stem.replace("ligand_", "")
    output_file = f"{output_prefix}_{ligand_name}_out.pdbqt"
    log_file = f"{output_prefix}_{ligand_name}_log.txt"
    
    print(f"リガンド {ligand_name} のドッキング計算を実行中...")
    
    # Vinaのコマンドを構築
    cmd = [
        VINA_EXECUTABLE,
        "--receptor", receptor_file,
        "--ligand", ligand_file,
        "--center_x", str(DOCKING_CENTER[0]),
        "--center_y", str(DOCKING_CENTER[1]),
        "--center_z", str(DOCKING_CENTER[2]),
        "--size_x", str(BOX_SIZE[0]),
        "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--out", output_file,
        "--log", log_file,
        "--exhaustiveness", str(EXHAUSTIVENESS),
        "--num_modes", str(NUM_MODES),
        "--cpu", str(CPU)
    ]
    
    try:
        # Vinaを実行
        start_time = time.time()
        result = subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        end_time = time.time()
        
        # 実行時間を計算
        elapsed_time = end_time - start_time
        
        print(f"ドッキング計算が完了しました: {output_file}")
        print(f"計算時間: {elapsed_time:.2f}秒")
        
        # ログファイルからベストスコアを抽出
        best_score = extract_best_score(log_file)
        print(f"ベストスコア: {best_score:.2f} kcal/mol")
        
        return output_file, best_score, log_file
    except subprocess.CalledProcessError as e:
        print(f"エラー: ドッキング計算に失敗しました: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None, None, None

def extract_best_score(log_file):
    """
    Vinaのログファイルからベストスコアを抽出
    
    Args:
        log_file: Vinaのログファイルパス
        
    Returns:
        ベストスコア（浮動小数点数）
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
            # モード行の3行後（テーブル最初のデータ行）を取得
            score_line = lines[mode_line_idx + 3]
            parts = score_line.strip().split()
            if len(parts) >= 2:
                # 2番目の列（affinity値）を取得
                return float(parts[1])
        
        return 0.0
    except Exception as e:
        print(f"警告: スコアの抽出に失敗しました: {e}")
        return 0.0

def plot_docking_results(results):
    """
    ドッキング結果をプロット
    
    Args:
        results: (リガンド名, スコア)のリスト
    """
    if not results:
        print("有効な結果がありません")
        return
    
    # スコアでソート
    sorted_results = sorted(results, key=lambda x: x[1])
    
    # プロット用にデータを準備
    ligand_names = [item[0] for item in sorted_results]
    scores = [item[1] for item in sorted_results]
    
    # プロット
    plt.figure(figsize=(10, 6))
    plt.barh(ligand_names, scores)
    plt.xlabel('結合エネルギー (kcal/mol)')
    plt.ylabel('リガンド')
    plt.title('ドッキング結果: 結合エネルギー比較')
    plt.tight_layout()
    
    # 保存
    plot_path = os.path.join(RESULTS_DIR, 'docking_results.png')
    plt.savefig(plot_path)
    print(f"プロットを保存しました: {plot_path}")
    
    return plot_path

def main():
    """メイン関数"""
    print(f"=== AutoDock Vinaによるドッキング計算 ===")
    print(f"レセプター: {RECEPTOR_PDBQT}")
    print(f"ドッキング中心: {DOCKING_CENTER}")
    print(f"ボックスサイズ: {BOX_SIZE}")
    print(f"Vina実行ファイル: {VINA_EXECUTABLE}")
    
    # レセプターファイルの確認
    if not os.path.exists(RECEPTOR_PDBQT):
        print(f"エラー: レセプターファイル {RECEPTOR_PDBQT} が見つかりません")
        sys.exit(1)
    
    # リガンドファイルの検索
    ligand_files = sorted(glob.glob(os.path.join(LIGANDS_DIR, "ligand_*.pdbqt")))
    if not ligand_files:
        print(f"エラー: {LIGANDS_DIR} にリガンドファイルが見つかりません")
        sys.exit(1)
    
    print(f"検出されたリガンド: {len(ligand_files)}ファイル")
    
    # 処理するリガンド数の制限
    ligand_files = ligand_files[:min(MAX_LIGANDS, len(ligand_files))]
    print(f"処理するリガンド: {len(ligand_files)}ファイル")
    
    # 結果を保存するリスト
    results = []
    
    # カレントディレクトリを保存
    original_cwd = os.getcwd()
    
    # 作業ディレクトリを結果ディレクトリに変更
    os.chdir(RESULTS_DIR)
    
    # 各リガンドに対してドッキングを実行
    for i, ligand_file in enumerate(ligand_files):
        ligand_name = Path(ligand_file).stem.replace("ligand_", "")
        print(f"\nリガンド {i+1}/{len(ligand_files)}: {ligand_name}")
        
        # 相対パスに変換
        rel_receptor = os.path.relpath(os.path.join(original_cwd, RECEPTOR_PDBQT))
        rel_ligand = os.path.relpath(os.path.join(original_cwd, ligand_file))
        
        # ドッキング計算
        output_file, score, log_file = run_vina_docking(
            rel_receptor, 
            rel_ligand, 
            "docking"
        )
        
        if output_file and score:
            results.append((ligand_name, score))
    
    # カレントディレクトリを元に戻す
    os.chdir(original_cwd)
    
    # 結果の表示
    print("\n===== ドッキング結果のランキング =====")
    print(f"{'順位':<5}{'リガンド':<20}{'スコア (kcal/mol)':<15}")
    print("-" * 40)
    
    sorted_results = sorted(results, key=lambda x: x[1])
    for i, (name, score) in enumerate(sorted_results):
        print(f"{i+1:<5}{name:<20}{score:<15.2f}")
    
    # 結果のプロット
    plot_path = plot_docking_results(results)
    
    print("\nすべてのドッキング計算が完了しました")
    print(f"結果ディレクトリ: {RESULTS_DIR}")

if __name__ == "__main__":
    main()