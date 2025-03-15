#!/usr/bin/env python3
"""
prepare_receptor コマンドを使用してPDBファイルをPDBQTに変換するスクリプト
"""

import os
import sys
import argparse
import logging
import subprocess
import tempfile
from pathlib import Path

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("prepare_protein")

def prepare_receptor(
    input_file: str,
    output_file: str = None,
    add_hydrogens: bool = True,
    preserve_charges: bool = False,
    cleanup: str = "nphs_lps_waters_nonstdres",
    remove_nonstd: bool = False,
    verbose: bool = False
) -> bool:
    """
    prepare_receptor コマンドを使用してタンパク質を準備する
    
    Args:
        input_file: 入力ファイルのパス（PDBファイル）
        output_file: 出力ファイルのパス（PDBQTファイル）
        add_hydrogens: 水素原子を追加するかどうか
        preserve_charges: 既存の電荷を保持するかどうか
        cleanup: クリーンアップの種類
        remove_nonstd: 標準でないアミノ酸残基を削除するかどうか
        verbose: 詳細な出力を表示するかどうか
        
    Returns:
        処理が成功したかどうか
    """
    try:
        # パスのチェック
        if not os.path.exists(input_file):
            logger.error(f"Input file does not exist: {input_file}")
            return False
        
        # 出力ファイルが指定されていない場合は入力ファイル名に基づいて生成
        if not output_file:
            input_path = Path(input_file)
            output_file = str(input_path.with_suffix('.pdbqt'))
        
        # prepare_receptor コマンドを構築
        cmd = [
            "prepare_receptor",
            "-r", input_file,
            "-o", output_file
        ]
        
        # オプションを追加
        if add_hydrogens:
            cmd.extend(["-A", "hydrogens"])
        
        if preserve_charges:
            cmd.append("-C")
        
        if cleanup:
            cmd.extend(["-U", cleanup])
        
        if remove_nonstd:
            cmd.extend(["-e", "True"])
        
        if verbose:
            cmd.append("-v")
        
        # コマンドを実行
        logger.info(f"Running prepare_receptor command: {' '.join(cmd)}")
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            check=False
        )
        
        # 結果をログに出力
        if result.stdout and verbose:
            logger.info(f"prepare_receptor stdout: {result.stdout}")
        
        if result.stderr:
            if result.returncode != 0:
                logger.error(f"prepare_receptor stderr: {result.stderr}")
            elif verbose:
                logger.warning(f"prepare_receptor stderr: {result.stderr}")
        
        # 結果を確認
        if result.returncode != 0:
            logger.error(f"prepare_receptor failed with return code {result.returncode}")
            return False
        
        # 出力ファイルが存在するか確認
        if not os.path.exists(output_file):
            logger.error(f"Output file was not created: {output_file}")
            return False
        
        logger.info(f"Protein preparation successful: {output_file}")
        return True
    
    except Exception as e:
        logger.error(f"Error in prepare_receptor: {e}")
        return False

def main():
    parser = argparse.ArgumentParser(description="Prepare protein using prepare_receptor")
    parser.add_argument("-r", "--receptor", required=True, help="Input receptor PDB file")
    parser.add_argument("-o", "--output", help="Output PDBQT file")
    parser.add_argument("-H", "--no-hydrogens", action="store_true", help="Do not add hydrogens")
    parser.add_argument("-C", "--preserve-charges", action="store_true", help="Preserve input charges")
    parser.add_argument("-U", "--cleanup", default="nphs_lps_waters_nonstdres", 
                        help="Cleanup types (nphs,lps,waters,nonstdres,deleteAltB)")
    parser.add_argument("-e", "--remove-nonstd", action="store_true", help="Remove non-standard residues")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose output")
    
    args = parser.parse_args()
    
    # 関数を呼び出す
    success = prepare_receptor(
        input_file=args.receptor,
        output_file=args.output,
        add_hydrogens=not args.no_hydrogens,
        preserve_charges=args.preserve_charges,
        cleanup=args.cleanup,
        remove_nonstd=args.remove_nonstd,
        verbose=args.verbose
    )
    
    # 終了コードを設定
    return 0 if success else 1

if __name__ == "__main__":
    sys.exit(main())