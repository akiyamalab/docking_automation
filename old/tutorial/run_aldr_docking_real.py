#!/usr/bin/env python
# coding: utf-8

"""
MeekoとVinaを使用してALDRタンパク質のドッキング計算を実行するスクリプト

このスクリプトは以下の処理を行います：
1. MeekoでタンパクとリガンドをPDBQT形式に変換
2. AutoDock Vinaでドッキング計算を実行
3. 結果を解析して保存
"""

import os
import sys
import time
import tempfile
import subprocess
import gzip
from pathlib import Path
from typing import List, Dict, Tuple, Any, Optional

import numpy as np
from rdkit import Chem
from meeko import MoleculePreparation
from meeko import PDBQTMolecule
from meeko import RDKitMolCreate

# ドッキング設定
DOCKING_CENTER = (16.611, -7.032, 14.448)  # 結晶リガンドから計算した中心座標
BOX_SIZE = (20.0, 20.0, 20.0)  # ボックスサイズ (x, y, z)
VINA_EXECUTABLE = "/advina/1.1.2/bin/vina"  # Vina実行ファイルのパス
MAX_COMPOUNDS = 5  # 処理する最大化合物数

# 作業ディレクトリの設定
WORK_DIR = Path(".")
INPUT_DIR = WORK_DIR / "input" / "ALDR"
OUTPUT_DIR = WORK_DIR / "output" / "ALDR"
os.makedirs(OUTPUT_DIR, exist_ok=True)

def prepare_receptor(pdb_file: Path) -> Path:
    """
    タンパク質のPDBファイルをPDBQT形式に変換
    
    Args:
        pdb_file: PDBファイルのパス
        
    Returns:
        PDBQT形式に変換されたファイルのパス
    """
    output_file = OUTPUT_DIR / f"{pdb_file.stem}.pdbqt"
    
    print(f"タンパク質 {pdb_file} をPDBQT形式に変換中...")
    
    # MeekoのPDBQTMoleculeを使用してPDBからPDBQTに変換
    # ただし、タンパク質の場合はMGLToolsのprepare_receptor.pyやOpen Babelを使うことが多い
    # ここでは簡易的にMeekoとmgltoolsのpython3実装を組み合わせて使用
    cmd = [
        "python3", "-m", "meeko.scripts.prepare_receptor",
        "-r", str(pdb_file),
        "-o", str(output_file)
    ]
    
    try:
        result = subprocess.run(
            cmd, 
            check=True, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"タンパク質の変換が完了しました: {output_file}")
        return output_file
    except subprocess.CalledProcessError as e:
        print(f"エラー: タンパク質の変換に失敗しました: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        
        # MGLToolsの代わりにOpen Babelを使ってみる
        print("Open Babelを使用してPDBQT変換を試みます...")
        cmd = [
            "obabel", 
            str(pdb_file), 
            "-O", str(output_file),
            "-xr"  # -xrオプションはPDBQT形式の変換のためのもの
        ]
        
        try:
            subprocess.run(cmd, check=True)
            print(f"Open Babelによるタンパク質の変換が完了しました: {output_file}")
            return output_file
        except subprocess.CalledProcessError:
            # どちらも失敗した場合、オプションとして直接ADTのscriptsフォルダのprepare_receptor4.pyを呼び出すこともできる
            # ただし、tutorial用なので簡易的な方法として、タンパク質の変換ができない場合はVinaの--receptorオプションで
            # PDBファイルを直接指定することもできる（非推奨だが可能）
            print(f"警告: タンパク質の変換に失敗しました。PDBファイルをそのまま使用します。")
            return pdb_file

def prepare_ligand(mol: Chem.Mol, name: str) -> Tuple[str, Optional[Path]]:
    """
    RDKitの分子オブジェクトをPDBQT形式に変換
    
    Args:
        mol: RDKitの分子オブジェクト
        name: リガンド名
        
    Returns:
        リガンド名とPDBQT形式のファイルパスのタプル
    """
    output_file = OUTPUT_DIR / f"ligand_{name}.pdbqt"
    temp_sdf = OUTPUT_DIR / f"temp_{name}.sdf"
    
    print(f"リガンド {name} をPDBQT形式に変換中...")
    
    try:
        # 一時SDFファイルに書き出し
        writer = Chem.SDWriter(str(temp_sdf))
        writer.write(mol)
        writer.close()
        
        # Open Babelを使用してPDBQT形式に変換
        cmd = [
            "obabel",
            str(temp_sdf),
            "-O", str(output_file),
            "-xh",  # 水素を追加
            "--partialcharge", "gasteiger",  # Gasteigerの部分電荷を計算
            "-xr"  # PDBQTとして出力
        ]
        
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        # 一時ファイルを削除
        if temp_sdf.exists():
            os.remove(temp_sdf)
        
        if output_file.exists() and output_file.stat().st_size > 0:
            print(f"リガンドの変換が完了しました: {output_file}")
            return name, output_file
        else:
            print(f"エラー: リガンド {name} のPDBQTファイルが生成されませんでした")
            return name, None
    except Exception as e:
        print(f"エラー: リガンド {name} の変換に失敗しました: {e}")
        # 一時ファイルを削除
        if temp_sdf.exists():
            os.remove(temp_sdf)
        return name, None

def run_vina_docking(receptor_file: Path, ligand_file: Path, ligand_name: str) -> Tuple[Optional[Path], Optional[float]]:
    """
    AutoDock Vinaを使用してドッキング計算を実行
    
    Args:
        receptor_file: レセプター（タンパク質）のPDBQTファイルパス
        ligand_file: リガンド（化合物）のPDBQTファイルパス
        ligand_name: リガンド名
        
    Returns:
        結果ファイルのパスとベストスコアのタプル
    """
    output_file = OUTPUT_DIR / f"docking_result_{ligand_name}.pdbqt"
    log_file = OUTPUT_DIR / f"docking_log_{ligand_name}.txt"
    
    print(f"リガンド {ligand_name} のドッキング計算を実行中...")
    
    # Vinaのコマンドを構築
    cmd = [
        VINA_EXECUTABLE,
        "--receptor", str(receptor_file),
        "--ligand", str(ligand_file),
        "--center_x", str(DOCKING_CENTER[0]),
        "--center_y", str(DOCKING_CENTER[1]),
        "--center_z", str(DOCKING_CENTER[2]),
        "--size_x", str(BOX_SIZE[0]),
        "--size_y", str(BOX_SIZE[1]),
        "--size_z", str(BOX_SIZE[2]),
        "--out", str(output_file),
        "--log", str(log_file),
        "--exhaustiveness", "8",  # 探索の徹底度
        "--num_modes", "9",       # 出力するポーズ数
        "--cpu", "1"              # 使用するCPU数
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
        
        return output_file, best_score
    except subprocess.CalledProcessError as e:
        print(f"エラー: ドッキング計算に失敗しました: {e}")
        print(f"STDOUT: {e.stdout}")
        print(f"STDERR: {e.stderr}")
        return None, None

def extract_best_score(log_file: Path) -> float:
    """
    Vinaのログファイルからベストスコアを抽出
    
    Args:
        log_file: Vinaのログファイルパス
        
    Returns:
        ベストスコア（浮動小数点数）
    """
    try:
        with open(log_file, 'r') as f:
            for line in f:
                if "mode |   affinity" in line:
                    # スコア行を探す
                    next_line = next(f, None)
                    if next_line:
                        # 最初のスコア（ベストスコア）を抽出
                        parts = next_line.split()
                        if len(parts) >= 3:
                            return float(parts[1])
        return 0.0
    except Exception as e:
        print(f"警告: スコアの抽出に失敗しました: {e}")
        return 0.0

def analyze_results(results: List[Tuple[str, float]]) -> None:
    """
    ドッキング結果を解析して表示
    
    Args:
        results: (リガンド名, スコア)のリスト
    """
    if not results:
        print("有効な結果がありません")
        return
    
    # スコアでソート
    sorted_results = sorted(results, key=lambda x: x[1])
    
    print("\n===== ドッキング結果のランキング =====")
    print(f"{'順位':<5}{'リガンド':<20}{'スコア (kcal/mol)':<15}")
    print("-" * 40)
    
    for i, (name, score) in enumerate(sorted_results):
        print(f"{i+1:<5}{name:<20}{score:<15.2f}")

def main():
    """メイン関数"""
    print(f"=== ALDRタンパク質とリガンドのドッキング計算 ===")
    print(f"ドッキング中心: {DOCKING_CENTER}")
    print(f"ボックスサイズ: {BOX_SIZE}")
    print(f"Vina実行ファイル: {VINA_EXECUTABLE}")
    print(f"入力ディレクトリ: {INPUT_DIR}")
    print(f"出力ディレクトリ: {OUTPUT_DIR}")
    
    # 入力ファイルのパス
    receptor_path = INPUT_DIR / "receptor.pdb"
    ligand_path = INPUT_DIR / "actives_final.sdf.gz"
    
    # タンパク質の準備
    receptor_pdbqt = prepare_receptor(receptor_path)
    
    # 圧縮されたSDFファイルを展開
    temp_sdf = OUTPUT_DIR / "actives_final.sdf"
    with gzip.open(ligand_path, 'rb') as f_in:
        with open(temp_sdf, 'wb') as f_out:
            f_out.write(f_in.read())
    
    # SDFファイルの読み込み
    suppl = Chem.SDMolSupplier(str(temp_sdf))
    mols = []
    for mol in suppl:
        if mol is not None:
            # 明示的な水素原子を追加（Meekoの要件）
            mol = Chem.AddHs(mol)
            mols.append(mol)
    print(f"読み込んだ化合物数: {len(mols)}")
    
    # 処理する化合物数の制限
    mols = mols[:min(MAX_COMPOUNDS, len(mols))]
    
    # 結果を保存するリスト
    results = []
    
    # 各化合物に対してドッキングを実行
    for i, mol in enumerate(mols):
        mol_name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"Compound_{i+1}"
        print(f"\nリガンド {i+1}/{len(mols)}: {mol_name} を処理中...")
        
        # リガンドの準備
        ligand_name, ligand_pdbqt = prepare_ligand(mol, mol_name)
        
        if ligand_pdbqt is None:
            print(f"リガンド {mol_name} の準備に失敗したため、スキップします")
            continue
        
        # ドッキング計算
        result_file, score = run_vina_docking(receptor_pdbqt, ligand_pdbqt, ligand_name)
        
        if result_file is not None and score is not None:
            results.append((ligand_name, score))
    
    # 結果の解析
    analyze_results(results)
    print("\nすべてのドッキング計算が完了しました")

if __name__ == "__main__":
    main()