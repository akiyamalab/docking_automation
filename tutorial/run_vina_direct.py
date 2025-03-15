#!/usr/bin/env python3
"""
既存のPDBQTファイルを使用して直接AutoDock Vinaでドッキング計算を実行するスクリプト
"""

import os
import sys
import time
import logging
import tempfile
import subprocess
import re
import shutil

# ロギング設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("vina_direct")

# 現在の作業ディレクトリを取得
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def run_vina_docking(ligand_path, receptor_path, output_dir,
                     center_x=0.0, center_y=0.0, center_z=0.0,
                     size_x=20.0, size_y=20.0, size_z=20.0,
                     exhaustiveness=8, num_modes=9):
    """AutoDock Vinaを使ってドッキング計算を実行"""
    os.makedirs(output_dir, exist_ok=True)
    
    # 出力ファイルのパス
    output_path = os.path.join(output_dir, "docking_result.pdbqt")
    log_path = os.path.join(output_dir, "docking_log.txt")
    
    # 一時設定ファイルの作成
    fd, config_path = tempfile.mkstemp(suffix=".conf", prefix="vina_")
    with os.fdopen(fd, 'w') as f:
        # グリッドボックスの設定
        f.write(f"center_x = {center_x}\n")
        f.write(f"center_y = {center_y}\n")
        f.write(f"center_z = {center_z}\n")
        f.write(f"size_x = {size_x}\n")
        f.write(f"size_y = {size_y}\n")
        f.write(f"size_z = {size_z}\n")
        
        # その他のパラメータ
        f.write(f"exhaustiveness = {exhaustiveness}\n")
        f.write(f"num_modes = {num_modes}\n")
        f.write("energy_range = 3\n")
    
    # Vinaパスの設定
    vina_path = "/advina/1.1.2/bin/vina"
    
    # Vinaコマンドの構築
    cmd = [
        vina_path,
        "--receptor", receptor_path,
        "--ligand", ligand_path,
        "--config", config_path,
        "--out", output_path,
        "--log", log_path,
        "--cpu", "1"
    ]
    
    try:
        # コマンドを実行
        logger.info(f"Running Vina command: {' '.join(cmd)}")
        start_time = time.time()
        
        process = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        
        end_time = time.time()
        execution_time = end_time - start_time
        
        # 結果を確認
        if process.returncode != 0:
            logger.error(f"Vina execution failed: {process.stderr}")
            return None
        
        logger.info(f"Vina calculation completed in {execution_time:.2f} seconds")
        
        # 結果を解析
        scores = parse_vina_log(log_path)
        
        return {
            "output_path": output_path,
            "log_path": log_path,
            "execution_time": execution_time,
            "scores": scores,
            "stdout": process.stdout,
            "stderr": process.stderr
        }
    
    except Exception as e:
        logger.error(f"Error executing Vina: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return None
    
    finally:
        # 一時ファイルの削除
        if os.path.exists(config_path):
            os.unlink(config_path)

def parse_vina_log(log_path):
    """Vinaのログファイルを解析してスコアを抽出"""
    scores = []
    
    try:
        with open(log_path, 'r') as f:
            log_content = f.read()
            
            # スコアテーブルの位置を特定
            table_start = log_content.find("-----+------------+----------+----------")
            if table_start == -1:
                logger.warning("Score table not found in log file")
                return scores
            
            # モードのスコアを抽出
            # 例: "   1       -8.7      0.000      0.000"
            mode_pattern = r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)"
            lines = log_content[table_start:].split('\n')
            
            for line in lines[1:]:  # ヘッダー行をスキップ
                match = re.match(mode_pattern, line.strip())
                if match:
                    mode_num = int(match.group(1))
                    affinity = float(match.group(2))
                    rmsd_lb = float(match.group(3))
                    rmsd_ub = float(match.group(4))
                    
                    scores.append({
                        'mode': mode_num,
                        'affinity': affinity,
                        'rmsd_lb': rmsd_lb,
                        'rmsd_ub': rmsd_ub
                    })
    
    except Exception as e:
        logger.error(f"Error parsing log file: {e}")
    
    return scores

def main():
    # 入力ファイルの確認
    ligand_path = os.path.join(CURRENT_DIR, "ligand.pdbqt")
    receptor_path = os.path.join(CURRENT_DIR, "protein.pdbqt")
    
    if not os.path.exists(ligand_path):
        logger.error(f"Ligand file not found: {ligand_path}")
        return 1
    
    if not os.path.exists(receptor_path):
        logger.error(f"Receptor file not found: {receptor_path}")
        return 1
    
    # 出力ディレクトリの作成
    output_dir = os.path.join(CURRENT_DIR, "direct_docking_results")
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    os.makedirs(output_dir, exist_ok=True)
    
    # グリッドボックスの設定
    # 通常はタンパク質構造の活性部位に合わせて設定します
    grid_center = (0.0, 0.0, 0.0)  # アプリケーションに応じて変更
    grid_size = (25.0, 25.0, 25.0)  # 一般的なサイズ
    
    # ドッキング計算を実行
    logger.info("Running docking calculation...")
    result = run_vina_docking(
        ligand_path=ligand_path,
        receptor_path=receptor_path,
        output_dir=output_dir,
        center_x=grid_center[0],
        center_y=grid_center[1],
        center_z=grid_center[2],
        size_x=grid_size[0],
        size_y=grid_size[1],
        size_z=grid_size[2],
        exhaustiveness=8,
        num_modes=9
    )
    
    if not result:
        logger.error("Docking calculation failed")
        return 1
    
    # 結果の要約を保存
    summary_path = os.path.join(output_dir, "docking_summary.txt")
    with open(summary_path, "w") as f:
        f.write("Docking Results Summary\n")
        f.write("======================\n\n")
        f.write(f"Ligand: {os.path.basename(ligand_path)}\n")
        f.write(f"Receptor: {os.path.basename(receptor_path)}\n")
        f.write(f"Execution time: {result['execution_time']:.2f} seconds\n\n")
        
        f.write("Binding Affinities:\n")
        f.write("-----------------\n")
        
        for score in result["scores"]:
            f.write(f"Mode {score['mode']}: {score['affinity']:.2f} kcal/mol (RMSD: {score['rmsd_lb']:.2f})\n")
        
        f.write("\nCommand Output:\n")
        f.write("---------------\n")
        f.write(result.get("stdout", "No output") + "\n")
        
        if result.get("stderr"):
            f.write("\nError Output:\n")
            f.write("-------------\n")
            f.write(result["stderr"] + "\n")
    
    logger.info(f"Results summary saved to: {summary_path}")
    
    # 結果ファイルをコピー
    shutil.copy(result["output_path"], os.path.join(output_dir, "best_pose.pdbqt"))
    shutil.copy(result["log_path"], os.path.join(output_dir, "vina_log.txt"))
    
    logger.info(f"Docking completed successfully! Results are in: {output_dir}")
    logger.info(f"Best binding affinity: {result['scores'][0]['affinity']:.2f} kcal/mol")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())