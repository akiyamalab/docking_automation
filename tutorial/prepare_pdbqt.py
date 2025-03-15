#!/usr/bin/env python
# coding: utf-8

"""
PDBQTファイル変換ユーティリティ

このスクリプトは、タンパク質（PDB形式）とリガンド（MOL2/SDF形式）を
AutoDock Vinaが認識可能なPDBQT形式に変換します。
変換にはOpen Babelを使用します。
"""

import os
import sys
import subprocess
import tempfile
from pathlib import Path
import gzip
from rdkit import Chem

class PDBQTConverter:
    """PDBQTファイル変換クラス"""
    
    def __init__(self, output_dir="./output"):
        """
        コンストラクタ
        
        Args:
            output_dir: 出力ディレクトリ
        """
        self.output_dir = Path(output_dir)
        os.makedirs(self.output_dir, exist_ok=True)
    
    def prepare_receptor(self, pdb_file):
        """
        タンパク質をPDBQT形式に変換
        
        Args:
            pdb_file: PDBファイルパス
            
        Returns:
            変換後のPDBQTファイルパス
        """
        pdb_path = Path(pdb_file)
        output_file = self.output_dir / f"{pdb_path.stem}.pdbqt"
        
        print(f"タンパク質 {pdb_file} をPDBQT形式に変換中...")
        
        try:
            # Open Babelを使用してPDBQTに変換
            cmd = [
                "obabel", 
                str(pdb_file), 
                "-O", str(output_file),
                "-xr"  # PDBQT形式として出力
            ]
            
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            
            if output_file.exists() and output_file.stat().st_size > 0:
                print(f"タンパク質の変換が完了しました: {output_file}")
                return output_file
            else:
                print(f"エラー: PDBQT変換に失敗しました")
                return None
                
        except subprocess.CalledProcessError as e:
            print(f"エラー: タンパク質の変換に失敗しました: {e}")
            print(f"STDOUT: {e.stdout}")
            print(f"STDERR: {e.stderr}")
            return None
    
    def prepare_ligand(self, ligand_file, ligand_id=None, temp_conv=False):
        """
        リガンドをPDBQT形式に変換
        
        Args:
            ligand_file: リガンドファイルパス（SDF, MOL2など）
            ligand_id: リガンド識別子（省略時はファイル名を使用）
            temp_conv: 一時ファイルにSDFを抽出してから変換するかどうか
            
        Returns:
            変換後のPDBQTファイルパスのリスト
        """
        ligand_path = Path(ligand_file)
        
        if ligand_id is None:
            ligand_id = ligand_path.stem
        
        print(f"リガンド {ligand_id} をPDBQT形式に変換中...")
        
        # gzファイルの場合は一時ファイルに展開
        if ligand_path.suffix == '.gz' or temp_conv:
            temp_dir = tempfile.mkdtemp()
            temp_file = Path(temp_dir) / ligand_path.stem
            
            if ligand_path.suffix == '.gz':
                with gzip.open(ligand_path, 'rb') as f_in:
                    with open(temp_file, 'wb') as f_out:
                        f_out.write(f_in.read())
                input_file = temp_file
            else:
                input_file = ligand_file
            
            # SDF/MOL2ファイルからリガンドを一つずつ抽出して変換
            output_files = self._convert_ligands_from_file(input_file)
            
            # 一時ファイル削除
            if temp_file.exists():
                os.remove(temp_file)
            os.rmdir(temp_dir)
            
            return output_files
        else:
            # 単一リガンドを変換
            output_file = self.output_dir / f"ligand_{ligand_id}.pdbqt"
            
            try:
                # Open Babelを使用してPDBQTに変換
                cmd = [
                    "obabel", 
                    str(ligand_file), 
                    "-O", str(output_file),
                    "-xh",  # 水素を追加
                    "--partialcharge", "gasteiger",  # Gasteigerの部分電荷を計算
                    "-p", "7.4",  # pH 7.4で電荷を計算
                    "-xr"  # PDBQT形式として出力
                ]
                
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                # PDBQT形式の修正（ROOT / TORSDOFタグの追加）
                self._fix_pdbqt_format(output_file)
                
                if output_file.exists() and output_file.stat().st_size > 0:
                    print(f"リガンド {ligand_id} の変換が完了しました: {output_file}")
                    return [output_file]
                else:
                    print(f"エラー: {ligand_id} のPDBQT変換に失敗しました")
                    return []
                    
            except subprocess.CalledProcessError as e:
                print(f"エラー: リガンド {ligand_id} の変換に失敗しました: {e}")
                print(f"STDOUT: {e.stdout}")
                print(f"STDERR: {e.stderr}")
                return []
    
    def _convert_ligands_from_file(self, input_file):
        """
        SDF/MOL2ファイルから複数のリガンドを抽出して変換
        
        Args:
            input_file: 入力ファイルパス
            
        Returns:
            変換後のPDBQTファイルパスのリスト
        """
        output_files = []
        
        # ファイル拡張子に基づいて処理
        path = Path(input_file)
        if path.suffix.lower() == '.sdf':
            # RDKitを使用してSDFからリガンドを読み込み
            suppl = Chem.SDMolSupplier(str(path))
            for i, mol in enumerate(suppl):
                if mol is not None:
                    # リガンド名の取得
                    if mol.HasProp("_Name"):
                        ligand_name = mol.GetProp("_Name")
                    else:
                        ligand_name = f"mol_{i+1}"
                    
                    # 一時SDFファイルに書き出し
                    temp_file = self.output_dir / f"temp_{ligand_name}.sdf"
                    writer = Chem.SDWriter(str(temp_file))
                    writer.write(mol)
                    writer.close()
                    
                    # PDBQTに変換
                    output_file = self.output_dir / f"ligand_{ligand_name}.pdbqt"
                    
                    try:
                        # Open Babelを使用してPDBQTに変換
                        cmd = [
                            "obabel", 
                            str(temp_file), 
                            "-O", str(output_file),
                            "-xh",  # 水素を追加
                            "--partialcharge", "gasteiger",  # Gasteigerの部分電荷を計算
                            "-p", "7.4",  # pH 7.4で電荷を計算
                            "-xr"  # PDBQT形式として出力
                        ]
                        
                        subprocess.run(cmd, check=True, capture_output=True, text=True)
                        
                        # PDBQT形式の修正（ROOT / TORSDOFタグの追加）
                        self._fix_pdbqt_format(output_file)
                        
                        if output_file.exists() and output_file.stat().st_size > 0:
                            print(f"リガンド {ligand_name} の変換が完了しました: {output_file}")
                            output_files.append(output_file)
                        
                        # 一時ファイルを削除
                        if temp_file.exists():
                            os.remove(temp_file)
                            
                    except subprocess.CalledProcessError as e:
                        print(f"エラー: リガンド {ligand_name} の変換に失敗しました: {e}")
                        if temp_file.exists():
                            os.remove(temp_file)
        
        elif path.suffix.lower() == '.mol2':
            # Open Babelを使用してMOL2ファイルを分割
            try:
                # 一時ディレクトリを作成
                temp_dir = tempfile.mkdtemp()
                temp_path = Path(temp_dir)
                temp_prefix = temp_path / "ligand"
                
                # MOL2ファイルを分割
                cmd = [
                    "obabel", 
                    str(path), 
                    "-O", f"{temp_prefix}*.mol2",
                    "-m"  # 複数分子に分割
                ]
                
                subprocess.run(cmd, check=True, capture_output=True, text=True)
                
                # 分割されたファイルをPDBQTに変換
                for mol2_file in temp_path.glob("*.mol2"):
                    ligand_name = mol2_file.stem
                    output_file = self.output_dir / f"ligand_{ligand_name}.pdbqt"
                    
                    try:
                        # Open Babelを使用してPDBQTに変換
                        cmd = [
                            "obabel", 
                            str(mol2_file), 
                            "-O", str(output_file),
                            "-xh",  # 水素を追加
                            "--partialcharge", "gasteiger",  # Gasteigerの部分電荷を計算
                            "-p", "7.4",  # pH 7.4で電荷を計算
                            "-xr"  # PDBQT形式として出力
                        ]
                        
                        subprocess.run(cmd, check=True, capture_output=True, text=True)
                        
                        # PDBQT形式の修正（ROOT / TORSDOFタグの追加）
                        self._fix_pdbqt_format(output_file)
                        
                        if output_file.exists() and output_file.stat().st_size > 0:
                            print(f"リガンド {ligand_name} の変換が完了しました: {output_file}")
                            output_files.append(output_file)
                        
                        # 一時ファイルを削除
                        os.remove(mol2_file)
                            
                    except subprocess.CalledProcessError as e:
                        print(f"エラー: リガンド {ligand_name} の変換に失敗しました: {e}")
                
                # 一時ディレクトリを削除
                os.rmdir(temp_dir)
                    
            except subprocess.CalledProcessError as e:
                print(f"エラー: MOL2ファイルの分割に失敗しました: {e}")
        
        return output_files
    
    def _fix_pdbqt_format(self, pdbqt_file):
        """
        PDBQT形式を修正（Vinaが認識可能な形式に）
        
        Args:
            pdbqt_file: PDBQTファイルパス
        """
        try:
            with open(pdbqt_file, 'r') as f:
                lines = f.readlines()
            
            # バックアップの作成
            with open(f"{pdbqt_file}.bak", 'w') as f:
                f.writelines(lines)
            
            # ATOM/HETATMの行数をカウント
            atom_count = sum(1 for line in lines if line.startswith('ATOM') or line.startswith('HETATM'))
            
            # ROOT/TORSDOFタグがあるか確認
            has_root = any('ROOT' in line for line in lines)
            has_torsdof = any('TORSDOF' in line for line in lines)
            
            if not has_root or not has_torsdof:
                # ROOTとTORSDOFタグを追加
                with open(pdbqt_file, 'w') as f:
                    f.write("ROOT\n")
                    for line in lines:
                        f.write(line)
                    f.write("ENDROOT\n")
                    f.write(f"TORSDOF {max(1, atom_count // 10)}\n")  # 原子数に基づいて回転可能結合数を推定
            
        except Exception as e:
            print(f"エラー: PDBQT形式の修正に失敗しました: {e}")

def main():
    """メイン関数"""
    if len(sys.argv) < 2:
        print("使用法: python prepare_pdbqt.py <receptor.pdb> [<ligand.sdf>...]")
        sys.exit(1)
    
    receptor_file = sys.argv[1]
    ligand_files = sys.argv[2:] if len(sys.argv) > 2 else []
    
    output_dir = Path("./output/ALDR")
    converter = PDBQTConverter(str(output_dir))
    
    # タンパク質の変換
    receptor_pdbqt = converter.prepare_receptor(receptor_file)
    
    # リガンドの変換
    ligand_pdbqts = []
    for ligand_file in ligand_files:
        ligand_pdbqts.extend(converter.prepare_ligand(ligand_file, temp_conv=True))
    
    print("\n変換結果:")
    print(f"タンパク質PDBQT: {receptor_pdbqt}")
    print(f"リガンドPDBQT: {len(ligand_pdbqts)} ファイル")
    
    return receptor_pdbqt, ligand_pdbqts

if __name__ == "__main__":
    main()