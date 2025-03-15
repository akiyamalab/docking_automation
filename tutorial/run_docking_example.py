#!/usr/bin/env python3
"""
新しいインポート構造を使用したドッキング自動化パッケージの使用例

このスクリプトでは、簡略化されたインポートパスを使用して：
1. 分子とタンパク質を読み込み
2. ドッキング計算を準備
3. モック実行を行う
"""

import os
import sys
import logging
from typing import Dict, Any

# 簡略化された新しいインポート構造を使用
from docking_automation.molecule import (
    Compound, 
    Protein,
    FormatType,
    MoleculeFormat,
    MoleculeStructure,
    Atom,
    Bond,
    MoleculeProperty
)

from docking_automation.docking import (
    Ligand,
    Receptor,
    GridBox,
    DockingTask,
    DockingConfiguration,
    DockingParameters,
    DockingParameter,
    ParameterType,
    TaskStatus
)

# ロギング設定
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("docking_example")

# 作業ディレクトリ取得
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def create_mock_protein() -> Protein:
    """モックタンパク質をメモリ上に作成"""
    # 簡易的な原子リスト（実際には多数の原子が存在）
    atoms = [
        Atom(atom_id=1, element="N", x=0.0, y=0.0, z=0.0),
        Atom(atom_id=2, element="C", x=1.0, y=0.0, z=0.0),
        Atom(atom_id=3, element="C", x=1.0, y=1.0, z=0.0),
        Atom(atom_id=4, element="O", x=2.0, y=1.0, z=0.0)
    ]
    
    # 簡易的な結合リスト
    bonds = [
        Bond(atom1_id=1, atom2_id=2, bond_type="single"),
        Bond(atom1_id=2, atom2_id=3, bond_type="single"),
        Bond(atom1_id=3, atom2_id=4, bond_type="double")
    ]
    
    # 構造と分子を作成
    structure = MoleculeStructure(atoms=atoms, bonds=bonds)
    
    # 以前は：docking_automation.domain.molecule.entity.protein.Protein
    # 新しい構造では：docking_automation.molecule.Protein
    protein = Protein(
        id="mock_protein",
        structure=structure,
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/protein.pdbqt",  # 仮想パス
        metadata={"name": "サンプルタンパク質"},  # 名前をメタデータに設定
        is_prepared=True
    )
    
    return protein

def create_mock_compound() -> Compound:
    """モック化合物をメモリ上に作成"""
    # アスピリンに似た簡易的な原子リスト
    atoms = [
        Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
        Atom(atom_id=2, element="C", x=1.0, y=0.0, z=0.0),
        Atom(atom_id=3, element="O", x=1.0, y=1.0, z=0.0),
        Atom(atom_id=4, element="O", x=2.0, y=0.0, z=0.0)
    ]
    
    # 簡易的な結合リスト
    bonds = [
        Bond(atom1_id=1, atom2_id=2, bond_type="single"),
        Bond(atom1_id=2, atom2_id=3, bond_type="single"),
        Bond(atom1_id=2, atom2_id=4, bond_type="double")
    ]
    
    # 構造と分子を作成
    structure = MoleculeStructure(atoms=atoms, bonds=bonds)
    
    # 以前は：docking_automation.domain.molecule.entity.compound.Compound
    # 新しい構造では：docking_automation.molecule.Compound
    compound = Compound(
        id="aspirin",
        structure=structure,
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/aspirin.pdbqt",  # 仮想パス
        metadata={"name": "アスピリン"},  # 名前をメタデータに設定
        is_prepared=True
    )
    
    return compound

def setup_docking_task() -> DockingTask:
    """ドッキングタスクをセットアップ"""
    # タンパク質とレセプターの作成
    protein = create_mock_protein()
    
    # 以前は：docking_automation.domain.docking.entity.receptor.Receptor
    # 新しい構造では：docking_automation.docking.Receptor
    receptor = Receptor(protein=protein)
    
    # 化合物とリガンドの作成
    compound = create_mock_compound()
    
    # 以前は：docking_automation.domain.docking.entity.ligand.Ligand
    # 新しい構造では：docking_automation.docking.Ligand
    ligand = Ligand(compound=compound)
    
    # グリッドボックスの設定
    # 以前は：docking_automation.domain.docking.value_object.grid_box.GridBox
    # 新しい構造では：docking_automation.docking.GridBox
    grid_box = GridBox(
        center_x=0.0, 
        center_y=0.0, 
        center_z=0.0,
        size_x=20.0, 
        size_y=20.0, 
        size_z=20.0
    )
    
    # ドッキングパラメータの設定
    # 以前は：docking_automation.domain.docking.value_object.docking_parameter.DockingParameters
    # 新しい構造では：docking_automation.docking.DockingParameters
    parameters = DockingParameters()
    parameters.add(DockingParameter(name="exhaustiveness", value=8, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="num_modes", value=9, parameter_type=ParameterType.INTEGER))
    
    # ドッキング設定の作成
    # 以前は：docking_automation.domain.docking.value_object.docking_configuration.DockingConfiguration
    # 新しい構造では：docking_automation.docking.DockingConfiguration
    configuration = DockingConfiguration(
        grid_box=grid_box,
        parameters=parameters,
        name="Default Configuration"
    )
    
    # ドッキングタスクの作成
    # 以前は：docking_automation.domain.docking.entity.docking_task.DockingTask
    # 新しい構造では：docking_automation.docking.DockingTask
    task = DockingTask(
        id="task_123",
        ligand=ligand,
        receptor=receptor,
        configuration=configuration,
        status=TaskStatus.PENDING
    )
    
    # デモンストレーション用の説明用データ
    # DockingTaskのreceptorプロパティは実際にはReceptorオブジェクトですが
    # 型チェックの問題を避けるため、メタデータを直接設定します
    # （このコードは実行時には問題ありません）
    task.metadata = {
        "receptor_name": protein.metadata.get("name", protein.id),
        "ligand_name": compound.metadata.get("name", compound.id)
    }
    
    return task

def main():
    """メイン関数"""
    logger.info("新しいインポート構造を使用したドッキング自動化パッケージのデモンストレーション")
    
    try:
        # タスクをセットアップ
        task = setup_docking_task()
        logger.info(f"タスク '{task.id}' を作成しました")
        
        # リガンド名とレセプター名を取得（安全な方法で）
        ligand_name = task.metadata.get("ligand_name", "不明なリガンド")
        receptor_name = task.metadata.get("receptor_name", "不明なレセプター")
        
        logger.info(f"リガンド: {ligand_name}")
        logger.info(f"レセプター: {receptor_name}")
        
        # タスクの準備状態を確認
        if task.is_ready():
            logger.info("タスクは実行可能です")
            task.mark_as_running()
            logger.info(f"タスクの状態: {task.status.name}")
            
            # 実際のドッキング計算はここで行われます
            # （このデモでは実行しません）
            
            task.mark_as_completed()
            logger.info(f"タスクの状態: {task.status.name}")
        else:
            logger.warning("タスクは実行できません")
            if task.error_message:
                logger.warning(f"エラー: {task.error_message}")
        
        logger.info("デモンストレーション完了")
        return 0
        
    except Exception as e:
        logger.error(f"エラーが発生しました: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())