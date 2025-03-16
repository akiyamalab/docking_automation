#!/usr/bin/env python3
"""
Meekoのインターフェースを確認するための簡単なスクリプト
"""

import sys
from rdkit import Chem
from meeko import MoleculePreparation, PDBQTWriterLegacy

def test_meeko():
    # 単純な分子を作成
    mol = Chem.MolFromSmiles('CC(=O)OC1=CC=CC=C1C(=O)O')
    if mol is None:
        print("Failed to create molecule")
        return
    
    # 3D構造を生成
    mol = Chem.AddHs(mol)
    print("Molecule created with", mol.GetNumAtoms(), "atoms")
    
    # MoleculePreparationクラスのプロパティとメソッドを調査
    prep = MoleculePreparation()
    print("\nMoleculePreparation dir:", dir(prep))
    
    # 分子を準備
    prep.prepare(mol)
    
    # 準備後のプロパティを確認
    print("\nAfter prepare, available properties:")
    for attr in dir(prep):
        if not attr.startswith("_"):
            print(" -", attr)
    
    # PDBQTWriterLegacyクラスのプロパティとメソッドを調査
    writer = PDBQTWriterLegacy()
    print("\nPDBQTWriterLegacy dir:", dir(writer))
    
    # 一時ファイルにPDBQTを書き出し
    output_file = "test_output.pdbqt"
    try:
        # 使用可能なメソッドを確認して最適なものを使用
        if hasattr(writer, "write"):
            writer.write(prep, output_file)
            print(f"\nWrote PDBQT using writer.write(prep, {output_file})")
        elif hasattr(writer, "write_pdbqt_file"):
            try:
                writer.write_pdbqt_file(prep.mol_info, output_file)
                print(f"\nWrote PDBQT using writer.write_pdbqt_file(prep.mol_info, {output_file})")
            except AttributeError:
                if hasattr(prep, "setup"):
                    writer.write_pdbqt_file(prep.setup, output_file)
                    print(f"\nWrote PDBQT using writer.write_pdbqt_file(prep.setup, {output_file})")
                else:
                    print("\nCan't identify correct attribute. Available prep attributes:")
                    for attr in dir(prep):
                        if not attr.startswith("_"):
                            print(f" - {attr}: {type(getattr(prep, attr))}")
    except Exception as e:
        print(f"Error writing PDBQT file: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    test_meeko()