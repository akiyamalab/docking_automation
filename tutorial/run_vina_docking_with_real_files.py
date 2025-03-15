#!/usr/bin/env python3
"""
実際のVinaコマンドと実際のファイルを使用してドッキング計算を実行するスクリプト
"""

import os
import sys
import logging
import time
import uuid
from pathlib import Path

# 自作モジュールへのパスを追加
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("vina_docking")

# 必要なドメインクラスをインポート
from docking_automation.domain.molecule.entity.compound import Compound
from docking_automation.domain.molecule.entity.protein import Protein
from docking_automation.domain.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.domain.molecule.value_object.molecule_structure import MoleculeStructure
from docking_automation.domain.docking.entity.ligand import Ligand
from docking_automation.domain.docking.entity.receptor import Receptor
from docking_automation.domain.docking.value_object.grid_box import GridBox
from docking_automation.domain.docking.value_object.docking_parameter import DockingParameters, DockingParameter, ParameterType

# Vinaサービスをインポート
from docking_automation.infrastructure.tools.vina.vina_docking_service import VinaDockingService


def main():
    # 出力ディレクトリの作成
    output_dir = os.path.join(current_dir, "vina_output")
    os.makedirs(output_dir, exist_ok=True)
    
    # 入力ファイルの確認
    ligand_path = os.path.join(current_dir, "ligand.pdbqt")
    receptor_path = os.path.join(current_dir, "protein.pdbqt")
    
    if not os.path.exists(ligand_path):
        logger.error(f"Ligand file not found: {ligand_path}")
        return 1
    
    if not os.path.exists(receptor_path):
        logger.error(f"Receptor file not found: {receptor_path}")
        return 1
    
    # Vinaパスを設定
    vina_path = "/advina/1.1.2/bin/vina"
    logger.info(f"Using Vina at: {vina_path}")
    
    # サービスのインスタンス化
    docking_service = VinaDockingService(vina_path=vina_path)
    
    # 1. 化合物とタンパク質のエンティティ作成（すでに準備済みのファイルを使用）
    logger.info(f"Loading ligand: {ligand_path}")
    ligand_format = MoleculeFormat(type=FormatType.PDBQT)
    compound = Compound(
        id=f"ligand_{uuid.uuid4().hex[:8]}",
        structure=MoleculeStructure(atoms=[], bonds=[]),  # ダミー構造（実際にはファイルから読み込む）
        format=ligand_format,
        path=ligand_path,
        is_prepared=True  # 既に準備済み
    )
    
    logger.info(f"Loading receptor: {receptor_path}")
    receptor_format = MoleculeFormat(type=FormatType.PDBQT)
    protein = Protein(
        id=f"receptor_{uuid.uuid4().hex[:8]}",
        structure=MoleculeStructure(atoms=[], bonds=[]),  # ダミー構造（実際にはファイルから読み込む）
        format=receptor_format,
        path=receptor_path,
        chains={"A"},
        is_prepared=True  # 既に準備済み
    )
    
    # 2. リガンドとレセプターの作成
    ligand = Ligand(compound=compound)
    receptor = Receptor(protein=protein)
    
    # 3. ドッキング設定の作成
    # タンパク質の中心にグリッドボックスを配置
    grid_box = GridBox(
        center_x=0.0,
        center_y=0.0,
        center_z=0.0,
        size_x=20.0,
        size_y=20.0,
        size_z=20.0
    )
    
    # ドッキングパラメータ
    parameters = DockingParameters()
    parameters.add(DockingParameter(name="exhaustiveness", value=8, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="num_modes", value=9, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="energy_range", value=3, parameter_type=ParameterType.INTEGER))
    parameters.add(DockingParameter(name="cpu", value=1, parameter_type=ParameterType.INTEGER))
    
    # デフォルト設定をベースにグリッドボックスを更新
    configuration = docking_service.get_default_configuration().with_grid_box(grid_box)
    
    # 4. ドッキングタスクの作成
    logger.info("Creating docking task")
    task = docking_service.create_task(
        ligand=ligand,
        receptor=receptor,
        configuration=configuration,
        metadata={"output_dir": output_dir}
    )
    
    # 5. ドッキング計算の実行
    logger.info(f"Executing docking task: {task.id}")
    try:
        result = docking_service.execute(task)
        
        # 6. 結果の表示
        logger.info("Docking calculation completed")
        logger.info(f"Task ID: {result.task.id}")
        logger.info(f"Execution time: {result.execution_time:.2f} seconds")
        logger.info(f"Number of poses: {len(result.poses)}")
        
        best_score = result.get_best_score()
        if best_score:
            logger.info(f"Best score: {best_score.value} {best_score.unit if best_score.unit else ''}")
        
        # 7. 結果の保存
        result_path = os.path.join(output_dir, "docking_result.txt")
        with open(result_path, "w") as f:
            f.write(f"Docking Result Summary\n")
            f.write(f"======================\n\n")
            f.write(f"Task ID: {result.task.id}\n")
            f.write(f"Ligand: {result.task.ligand.name}\n")
            # 複数レセプター対応
            try:
                if hasattr(result.task, 'get_primary_receptor') and callable(getattr(result.task, 'get_primary_receptor')):
                    # 新しいインターフェースを使用
                    receptor = result.task.get_primary_receptor()
                    receptor_info = str(receptor.name) if hasattr(receptor, 'name') else "Unknown"
                    
                    if hasattr(result.task, 'is_multi_receptor') and callable(getattr(result.task, 'is_multi_receptor')) and result.task.is_multi_receptor():
                        receptors = result.task.get_receptors()
                        receptor_info = f"{receptor_info} (+{len(receptors)-1} more)"
                    
                    f.write(f"Receptor: {receptor_info}\n")
                else:
                    # 互換性のためのフォールバック
                    receptor = result.task.receptor
                    if isinstance(receptor, list):
                        if receptor:
                            receptor_name = str(receptor[0].name) if hasattr(receptor[0], 'name') else "Unknown"
                            if len(receptor) > 1:
                                receptor_name = f"{receptor_name} (+{len(receptor)-1} more)"
                        else:
                            receptor_name = "None"
                    else:
                        receptor_name = str(receptor.name) if hasattr(receptor, 'name') else "Unknown"
                    
                    f.write(f"Receptor: {receptor_name}\n")
            except Exception as e:
                logger.warning(f"Failed to get receptor info: {e}")
                f.write("Receptor: Unknown\n")
            f.write(f"Poses: {len(result.poses)}\n")
            
            if best_score:
                f.write(f"Best Score: {best_score.value:.2f} {best_score.unit if best_score.unit else ''}\n")
                
            f.write(f"Execution Time: {result.execution_time:.2f} seconds\n\n")
            
            # ポーズごとのスコア
            f.write(f"Poses and Scores:\n")
            f.write(f"-----------------\n")
            for i, pose in enumerate(result.poses):
                f.write(f"Pose {pose.rank}: ")
                
                # ポーズのスコア
                pose_scores = result.get_scores_for_pose(pose.rank)
                if pose_scores:
                    score_str = ", ".join([f"{s.name}: {s.value:.2f}" for s in pose_scores])
                    f.write(f"{score_str}")
                else:
                    # ポーズごとのスコアがない場合は一般的なスコアを使用
                    for score in result.scores:
                        if score.name == "Binding Affinity":
                            f.write(f"Binding Affinity: {score.value:.2f} {score.unit if score.unit else ''}")
                            break
                
                f.write(f", RMSD to best: {pose.rmsd_to_best:.2f}\n")
            
        logger.info(f"Results saved to: {result_path}")
        return 0
        
    except Exception as e:
        logger.error(f"Error executing docking task: {e}")
        import traceback
        logger.error(traceback.format_exc())
        return 1


if __name__ == "__main__":
    sys.exit(main())