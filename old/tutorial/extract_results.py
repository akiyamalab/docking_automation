#!/usr/bin/env python
# coding: utf-8

"""
ドッキング結果を集計して表示するスクリプト

既存のVinaのログファイルからスコアを抽出して結果を表示します。
"""

import os
import sys
import glob
from pathlib import Path
import matplotlib
# GUIが不要なバックエンドに設定（サーバー環境用）
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# 結果ディレクトリ
RESULTS_DIR = "output/ALDR/results"

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
        
        # スコア行を直接探す
        for i, line in enumerate(lines):
            if line.strip().startswith('1 '):
                parts = line.strip().split()
                if len(parts) >= 2:
                    # 2番目の要素がスコア
                    try:
                        return float(parts[1])
                    except ValueError:
                        pass
        
        # 第2の方法: テーブル構造を解析
        # "mode |   affinity"行を探す
        for i, line in enumerate(lines):
            if "mode |   affinity" in line:
                # 区切り線をスキップ
                if i + 3 < len(lines):
                    # 最初のデータ行を取得
                    score_line = lines[i + 3]
                    parts = score_line.strip().split()
                    if len(parts) >= 2:
                        return float(parts[1])
        
        print(f"ファイル {log_file} からスコアを抽出できませんでした")
        return 0.0
    except Exception as e:
        print(f"警告: スコアの抽出に失敗しました: {e}")
        return 0.0

def plot_docking_results(results, output_path):
    """
    ドッキング結果をプロット
    
    Args:
        results: (リガンド名, スコア)のリスト
        output_path: 出力ファイルパス
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
    plt.savefig(output_path)
    print(f"プロットを保存しました: {output_path}")
    
    return output_path

def main():
    """メイン関数"""
    print(f"=== ドッキング結果の集計 ===")
    
    # ログファイルの検索
    log_files = sorted(glob.glob(os.path.join(RESULTS_DIR, "docking_*_log.txt")))
    if not log_files:
        print(f"エラー: {RESULTS_DIR} にログファイルが見つかりません")
        sys.exit(1)
    
    print(f"検出されたログファイル: {len(log_files)}ファイル")
    
    # 結果を保存するリスト
    results = []
    
    # 各ログファイルからスコアを抽出
    for log_file in log_files:
        log_path = Path(log_file)
        ligand_name = log_path.stem.split('_')[1]  # docking_CHEMBL123_log -> CHEMBL123
        
        # スコアの抽出
        score = extract_best_score(log_file)
        print(f"リガンド {ligand_name}: スコア = {score:.2f} kcal/mol")
        
        if score != 0.0:
            results.append((ligand_name, score))
    
    # 結果の表示
    print("\n===== ドッキング結果のランキング =====")
    print(f"{'順位':<5}{'リガンド':<20}{'スコア (kcal/mol)':<15}")
    print("-" * 40)
    
    sorted_results = sorted(results, key=lambda x: x[1])
    for i, (name, score) in enumerate(sorted_results):
        print(f"{i+1:<5}{name:<20}{score:<15.2f}")
    
    # 結果のプロット
    if results:
        plot_path = os.path.join(RESULTS_DIR, "docking_results.png")
        plot_docking_results(results, plot_path)
    
    print("\n集計が完了しました")

if __name__ == "__main__":
    main()