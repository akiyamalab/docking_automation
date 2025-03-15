#!/usr/bin/env python3
"""
新しいインポート構造を使用したドッキングワークフローの例

このスクリプトでは実際のドッキングワークフローを簡略化したデモを行います：
1. 複数の分子とタンパク質を読み込み
2. 準備とドッキング計算をセットアップ
3. 結果を解析して表示
"""

import os
import sys
import logging
import time
import uuid
from typing import Dict, Any, List

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
    DockingResult,
    DockingConfiguration,
    DockingParameters,
    DockingParameter,
    ParameterType,
    TaskStatus,
    Score,
    ScoreType
)

# ロギング設定
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("docking_workflow")

# 作業ディレクトリ取得
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))

def create_mock_compounds() -> List[Compound]:
    """モック化合物のリストを作成"""
    compounds = []
    
    # アスピリン
    aspirin = Compound(
        id="aspirin",
        structure=MoleculeStructure(atoms=[
            Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="O", x=1.0, y=0.0, z=0.0)
        ], bonds=[
            Bond(atom1_id=1, atom2_id=2, bond_type="single")
        ]),
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/aspirin.pdbqt",
        metadata={"name": "アスピリン", "molecular_weight": 180.16},
        is_prepared=True
    )
    compounds.append(aspirin)
    
    # イブプロフェン
    ibuprofen = Compound(
        id="ibuprofen",
        structure=MoleculeStructure(atoms=[
            Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="C", x=1.0, y=0.0, z=0.0)
        ], bonds=[
            Bond(atom1_id=1, atom2_id=2, bond_type="single")
        ]),
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/ibuprofen.pdbqt",
        metadata={"name": "イブプロフェン", "molecular_weight": 206.29},
        is_prepared=True
    )
    compounds.append(ibuprofen)
    
    # パラセタモール
    paracetamol = Compound(
        id="paracetamol",
        structure=MoleculeStructure(atoms=[
            Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="N", x=1.0, y=0.0, z=0.0)
        ], bonds=[
            Bond(atom1_id=1, atom2_id=2, bond_type="single")
        ]),
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/paracetamol.pdbqt",
        metadata={"name": "パラセタモール", "molecular_weight": 151.16},
        is_prepared=True
    )
    compounds.append(paracetamol)
    
    return compounds

def create_mock_protein() -> Protein:
    """モックタンパク質を作成"""
    return Protein(
        id="protein_1",
        structure=MoleculeStructure(atoms=[
            Atom(atom_id=1, element="N", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="C", x=1.0, y=0.0, z=0.0)
        ], bonds=[
            Bond(atom1_id=1, atom2_id=2, bond_type="single")
        ]),
        format=MoleculeFormat(type=FormatType.PDBQT),
        path="/path/to/protein.pdbqt",
        metadata={"name": "サンプルタンパク質", "pdb_id": "1ABC"},
        is_prepared=True
    )

def create_mock_docking_result(task: DockingTask) -> DockingResult:
    """モックドッキング結果を作成"""
    # スコアの作成
    scores = [
        Score(
            value=-8.7,
            score_type=ScoreType.BINDING_AFFINITY,
            name="Binding Affinity",
            unit="kcal/mol"
        ),
        Score(
            value=-7.5,
            score_type=ScoreType.BINDING_FREE_ENERGY,
            name="Binding Free Energy",
            unit="kcal/mol"
        )
    ]
    
    # 結果の作成
    result = DockingResult(
        id=str(uuid.uuid4())[:8],  # UUIDを文字列に変換してからスライス
        task=task,
        poses=[],  # 簡略化のためポーズは空リスト
        scores=scores,
        created_at=time.time(),
        execution_time=2.5,  # 模擬実行時間（秒）
        tool_info={"tool": "Vina", "version": "1.2.3"}
    )
    
    return result

def simulate_docking_workflow():
    """ドッキングワークフローのシミュレーション"""
    logger.info("=== ドッキングワークフロー開始 ===")
    
    # 1. 化合物の読み込み
    logger.info("化合物を読み込み中...")
    compounds = create_mock_compounds()
    logger.info(f"{len(compounds)}個の化合物を読み込みました")
    
    # 2. タンパク質の読み込み
    logger.info("タンパク質を読み込み中...")
    protein = create_mock_protein()
    logger.info(f"タンパク質 '{protein.metadata.get('name')}' を読み込みました")
    
    # 3. レセプターの作成
    receptor = Receptor(protein=protein)
    
    # 4. グリッドボックスの設定
    grid_box = GridBox(
        center_x=0.0, center_y=0.0, center_z=0.0,
        size_x=20.0, size_y=20.0, size_z=20.0
    )
    
    # 5. ドッキングパラメータの設定
    parameters = DockingParameters()
    parameters.add(DockingParameter(name="exhaustiveness", value=8, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="num_modes", value=9, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="energy_range", value=3, parameter_type=ParameterType.FLOAT))
    
    # 6. ドッキング設定の作成
    configuration = DockingConfiguration(
        grid_box=grid_box,
        parameters=parameters,
        name="Standard Docking Configuration"
    )
    
    # 7. 各化合物に対するドッキング計算
    results = []
    
    for compound in compounds:
        compound_name = compound.metadata.get("name", compound.id)
        logger.info(f"\n*** 化合物 '{compound_name}' の処理開始 ***")
        
        # リガンドの作成
        ligand = Ligand(compound=compound)
        
        # ドッキングタスクの作成 - UUIDを文字列に変換してからスライス
        task = DockingTask(
            id=f"task_{compound.id}_{str(uuid.uuid4())[:6]}",
            ligand=ligand,
            receptor=receptor,
            configuration=configuration,
            status=TaskStatus.PENDING
        )
        
        logger.info(f"タスク '{task.id}' を作成しました")
        
        # タスクの実行準備
        if task.is_ready():
            logger.info("タスクは実行可能です")
            task.mark_as_running()
            logger.info(f"タスクの状態: {task.status.name}")
            
            # 実際のドッキング計算をシミュレート
            logger.info(f"ドッキング計算を実行中...")
            time.sleep(0.5)  # シミュレーション用の遅延
            
            # 結果の生成
            result = create_mock_docking_result(task)
            logger.info(f"ドッキング計算が完了しました")
            
            # タスクを完了に設定
            task.mark_as_completed()
            logger.info(f"タスクの状態: {task.status.name}")
            
            # 結果を表示
            best_score = result.get_best_score()
            if best_score:
                logger.info(f"最良スコア: {best_score.value:.2f} {best_score.unit}")
            
            results.append(result)
        else:
            logger.warning(f"タスク '{task.id}' は実行できません")
            if task.error_message:
                logger.warning(f"エラー: {task.error_message}")
    
    # 8. 結果の要約
    logger.info("\n=== ドッキング結果の要約 ===")
    logger.info(f"実行したタスク数: {len(results)}/{len(compounds)}")
    
    if results:
        logger.info("化合物のランキング（最良スコア順）:")
        
        # スコアでソート（値が低いほど良い）
        sorted_results = sorted(
            results, 
            key=lambda r: min([s.value for s in r.scores]) if r.scores else float('inf')
        )
        
        for i, result in enumerate(sorted_results, 1):
            compound_name = result.task.ligand.compound.metadata.get("name", result.task.ligand.compound.id)
            best_score = min([s.value for s in result.scores]) if result.scores else "N/A"
            
            if isinstance(best_score, float):
                logger.info(f"{i}. {compound_name}: {best_score:.2f} kcal/mol")
            else:
                logger.info(f"{i}. {compound_name}: {best_score}")
    
    logger.info("\n=== ドッキングワークフロー完了 ===")
    return results

def main():
    try:
        simulate_docking_workflow()
        return 0
    except Exception as e:
        logger.error(f"エラーが発生しました: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1

if __name__ == "__main__":
    sys.exit(main())