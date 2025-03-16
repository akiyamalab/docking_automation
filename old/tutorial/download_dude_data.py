#!/usr/bin/env python
# coding: utf-8

"""
DUD-E データベースから ALDR タンパク質とリガンド情報をダウンロードして保存するスクリプト
"""

import os
import sys
import gzip
import shutil
import urllib.request
from pathlib import Path

# 出力ディレクトリの設定
OUTPUT_DIR = Path("input/ALDR")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# DUD-E データベースの URL
BASE_URL = "http://dude.docking.org/targets/aldr"
RECEPTOR_URL = f"{BASE_URL}/receptor.pdb"
CRYSTAL_LIGAND_URL = f"{BASE_URL}/crystal_ligand.mol2"
ACTIVES_URL = f"{BASE_URL}/actives_final.sdf.gz"
DECOYS_URL = f"{BASE_URL}/decoys_final.sdf.gz"

def download_file(url, output_path):
    """
    指定したURLからファイルをダウンロード
    
    Args:
        url: ダウンロード元URL
        output_path: 保存先パス
    """
    print(f"ダウンロード中: {url} -> {output_path}")
    try:
        urllib.request.urlretrieve(url, output_path)
        print(f"ダウンロード完了: {output_path}")
        return True
    except Exception as e:
        print(f"エラー: ダウンロードに失敗しました: {e}")
        return False

def main():
    """メイン関数"""
    print("DUD-E から ALDR タンパク質とリガンド情報をダウンロードしています...")
    
    # レセプター（タンパク質）のダウンロード
    receptor_path = OUTPUT_DIR / "receptor.pdb"
    download_file(RECEPTOR_URL, receptor_path)
    
    # 結晶リガンドのダウンロード
    crystal_ligand_path = OUTPUT_DIR / "crystal_ligand.mol2"
    download_file(CRYSTAL_LIGAND_URL, crystal_ligand_path)
    
    # 活性化合物のダウンロード
    actives_path = OUTPUT_DIR / "actives_final.sdf.gz"
    download_file(ACTIVES_URL, actives_path)
    
    # デコイ化合物のダウンロード（オプション）
    decoys_path = OUTPUT_DIR / "decoys_final.sdf.gz"
    download_file(DECOYS_URL, decoys_path)
    
    print("\nすべてのファイルのダウンロードが完了しました。")
    print(f"保存先ディレクトリ: {OUTPUT_DIR}")
    
    # ダウンロードしたファイルの一覧を表示
    print("\nダウンロードしたファイル:")
    for file_path in OUTPUT_DIR.glob("*"):
        file_size = file_path.stat().st_size / 1024  # キロバイト単位
        print(f"- {file_path.name} ({file_size:.1f} KB)")

if __name__ == "__main__":
    main()