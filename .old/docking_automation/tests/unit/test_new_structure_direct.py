#!/usr/bin/env python3
"""新しいディレクトリ構造のテスト (直接インポート版)

このスクリプトは、domainディレクトリからmoleculeおよびdockingディレクトリへの
移行が適切に行われたことを確認するためのテストを実行します。
トップレベルのインポートではなく、直接のインポートを使用します。
"""

import os
import sys
import uuid
from pathlib import Path

# ライブラリのパスを追加
sys.path.insert(0, os.path.abspath('.'))

# 新しいパスの確認
print("1. 移行先ディレクトリの存在確認")
print(f"molecule/entity: {os.path.exists('docking_automation/molecule/entity')}")
print(f"molecule/service: {os.path.exists('docking_automation/molecule/service')}")

# Compoundクラスのテスト
print("\n2. Compoundクラスのテスト")
try:
    # 直接インポート
    from docking_automation.molecule.entity.compound import Compound
    
    # Compoundインスタンスの作成
    compound = Compound(
        id=str(uuid.uuid4()),
        path="input_ligands/aspirin.pdb"
    )
    
    print(f"Compound作成成功: {compound}")
    print(f"Compound名: {compound.name}")
    
    # メタデータ追加のテスト
    compound.add_metadata({"test_key": "test_value"})
    print(f"メタデータ追加: {compound.metadata}")
    
except Exception as e:
    print(f"Compoundクラスのテストに失敗: {e}")
    import traceback
    traceback.print_exc()

# Proteinクラスのテスト
print("\n3. Proteinクラスのテスト")
try:
    # 直接インポート
    from docking_automation.molecule.entity.protein import Protein
    
    # Proteinインスタンスの作成
    protein = Protein(
        id=str(uuid.uuid4()),
        path="protein.pdb",
        chains={"A", "B"}
    )
    
    print(f"Protein作成成功: {protein}")
    print(f"Protein名: {protein.name}")
    print(f"チェーン数: {len(protein.chains)}")
    
    # メタデータ追加のテスト
    protein.add_metadata({"active_site": {"x": 10, "y": 20, "z": 30}})
    print(f"アクティブサイト定義: {protein.has_active_site_defined()}")
    
except Exception as e:
    print(f"Proteinクラスのテストに失敗: {e}")
    import traceback
    traceback.print_exc()

# 分子準備サービスのテスト
print("\n4. 分子準備サービスのテスト")
try:
    # 直接インポート
    from docking_automation.molecule.service.molecule_preparation_service import MoleculePreparationService
    from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService
    
    # 直接インスタンス作成
    preparation_service = VinaPreparationService()
    print(f"サービス作成成功: {preparation_service.__class__.__name__}")
    
    # サポートされているメソッドを取得
    methods = preparation_service.get_supported_methods()
    print(f"サポートされているメソッド: {list(methods.keys())}")
    
    # 化合物の準備
    if 'compound' in locals():
        prepared_compound = preparation_service.prepare_ligand(compound)
        print(f"化合物の準備成功: {prepared_compound.metadata}")
    
    # タンパク質の準備
    if 'protein' in locals():
        prepared_protein = preparation_service.prepare_receptor(protein)
        print(f"タンパク質の準備成功: {prepared_protein.metadata}")
    
except Exception as e:
    print(f"分子準備サービスのテストに失敗: {e}")
    import traceback
    traceback.print_exc()

print("\nテスト完了。エラーがなければ、ディレクトリ構造の移行は成功しています。")