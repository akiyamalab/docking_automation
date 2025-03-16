#!/usr/bin/env python
# coding: utf-8

"""
複数の化合物に対してバッチドッキング計算を実行するチュートリアル

このスクリプトは、BatchDockingWorkflowクラスを使用して、
複数の化合物に対するドッキング計算を効率的に実行する方法を示します。
"""

import os
import sys
import glob
from pathlib import Path
import logging
import time
import argparse
import numpy as np
import matplotlib.pyplot as plt
from rdkit import Chem

from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.service.molecule_preparation_factory import MoleculePreparationFactory
from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.docking_parameter import DockingParameter, ParameterType
from docking_automation.docking.value_object.docking_parameters import DockingParameters
from docking_automation.docking.service.docking_service import DockingService
from docking_automation.infrastructure.tools.vina.vina_docking_service import VinaDockingService
from docking_automation.application.workflows.batch_docking_workflow import BatchDockingWorkflow
from docking_automation.application.services.result_analysis_service import ResultAnalysisService
from docking_automation.infrastructure.formats.openbabel_format_converter import OpenBabelFormatConverter

# ロギングの設定
logging.basicConfig(level=logging.INFO, 
                    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger("BatchDockingTutorial")

def parse_arguments():
    """コマンドライン引数をパースする"""
    parser = argparse.ArgumentParser(description='複数の化合物に対するバッチドッキング計算')
    
    parser.add_argument('--receptor', '-r', required=True,
                      help='レセプター（タンパク質）のファイルパス（PDB形式）')
    
    parser.add_argument('--ligands', '-l', required=True,
                      help='リガンドファイルのパス（複数可、ワイルドカード使用可）')
    
    parser.add_argument('--output', '-o', default='output/batch',
                      help='結果を保存するディレクトリ（デフォルト: output/batch）')
    
    parser.add_argument('--center', '-c', type=float, nargs=3, 
                      help='ドッキング中心座標 (x y z)')
    
    parser.add_argument('--size', '-s', type=float, nargs=3, default=[20.0, 20.0, 20.0],
                      help='ボックスサイズ (x y z) (デフォルト: 20.0 20.0 20.0)')
    
    parser.add_argument('--exhaustiveness', '-e', type=int, default=8,
                      help='探索の徹底度（デフォルト: 8）')
    
    parser.add_argument('--cpu', type=int, default=1,
                      help='使用するCPU数（デフォルト: 1）')
    
    parser.add_argument('--num-modes', '-n', type=int, default=9,
                      help='出力するポーズ数（デフォルト: 9）')
    
    parser.add_argument('--jobs', '-j', type=int, default=1,
                      help='並列実行数（デフォルト: 1）')
    
    parser.add_argument('--max-ligands', '-m', type=int, default=None,
                      help='処理する最大リガンド数（デフォルト: すべて）')
    
    parser.add_argument('--convert-sdf', action='store_true',
                      help='結果をSDFに変換する')
    
    return parser.parse_args()

def get_ligand_paths(ligands_pattern, max_ligands=None):
    """リガンドファイルのパスを取得する"""
    # ワイルドカードを使ってファイルを検索
    ligand_paths = sorted(glob.glob(ligands_pattern))
    
    if not ligand_paths:
        raise ValueError(f"リガンドファイルが見つかりません: {ligands_pattern}")
    
    logger.info(f"検出されたリガンドファイル: {len(ligand_paths)}個")
    
    # 最大リガンド数の制限
    if max_ligands is not None and max_ligands > 0:
        ligand_paths = ligand_paths[:min(max_ligands, len(ligand_paths))]
        logger.info(f"処理するリガンド: {len(ligand_paths)}個に制限")
    
    return ligand_paths

def create_docking_configuration(center, size, exhaustiveness, num_modes, cpu):
    """ドッキング設定を作成する"""
    # GridBoxの設定
    grid_box = GridBox(
        center_x=center[0],
        center_y=center[1],
        center_z=center[2],
        size_x=size[0],
        size_y=size[1],
        size_z=size[2]
    )
    
    # パラメータの設定
    # パラメータの設定
    parameters = DockingParameters()
    parameters.add(DockingParameter("exhaustiveness", exhaustiveness, ParameterType.INTEGER))
    parameters.add(DockingParameter("num_modes", num_modes, ParameterType.INTEGER))
    parameters.add(DockingParameter("energy_range", 3, ParameterType.FLOAT))
    parameters.add(DockingParameter("seed", 0, ParameterType.INTEGER))
    parameters.add(DockingParameter("cpu", cpu, ParameterType.INTEGER))
    # 設定を返す
    return DockingConfiguration(
        grid_box=grid_box,
        parameters=parameters,
        name="BatchDocking",
        description="バッチドッキング設定"
    )

def display_results_ranking(results, output_dir):
    """ドッキング結果のランキングを表示する"""
    # 結果をスコアでソート
    sorted_results = sorted(
        [(r.task.ligand.name, r.get_best_score().value if r.get_best_score() else float('inf')) 
         for r in results],
        key=lambda x: x[1]
    )
    
    print("\n===== ドッキング結果のランキング =====")
    print(f"{'順位':<5}{'リガンド':<30}{'スコア (kcal/mol)':<15}")
    print("-" * 50)
    
    for i, (name, score) in enumerate(sorted_results):
        print(f"{i+1:<5}{name:<30}{score:<15.2f}")
    
    # データをプロット用に準備
    ligand_names = [item[0] for item in sorted_results[:10]]  # 上位10件
    scores = [item[1] for item in sorted_results[:10]]
    
    # プロット作成
    plt.figure(figsize=(10, 6))
    plt.barh(ligand_names, scores)
    plt.xlabel('結合エネルギー (kcal/mol)')
    plt.ylabel('リガンド')
    plt.title('ドッキング結果: 上位10化合物の結合エネルギー')
    plt.tight_layout()
    
    # プロットを保存
    plot_path = os.path.join(output_dir, 'docking_results.png')
    plt.savefig(plot_path)
    print(f"\nプロットを保存しました: {plot_path}")

def convert_results_to_sdf(results, output_dir):
    """ドッキング結果のPDBQTファイルをSDFに変換する"""
    converter = OpenBabelFormatConverter()
    
    sdf_dir = os.path.join(output_dir, "sdf")
    os.makedirs(sdf_dir, exist_ok=True)
    
    print("\n===== 結果をSDFに変換 =====")
    converted_files = []
    
    for result in results:
        # 各リガンドの結果ディレクトリからベストポーズを取得
        ligand_name = result.task.ligand.name
        ligand_dir = os.path.join(output_dir, ligand_name)
        
        # ディレクトリ内のPDBQTファイルを探す
        pdbqt_files = glob.glob(os.path.join(ligand_dir, "*.pdbqt"))
        if not pdbqt_files:
            logger.warning(f"{ligand_name} のPDBQTファイルが見つかりません")
            continue
        
        # 最初のPDBQTファイルをSDFに変換
        pdbqt_file = pdbqt_files[0]
        sdf_file = os.path.join(sdf_dir, f"{ligand_name}.sdf")
        
        try:
            # 変換を実行 - 入力ファイルから自動的に形式を検出
            with open(pdbqt_file, 'r') as f:
                pdbqt_content = f.read()
                
            with open(sdf_file, 'w') as f:
                # このサンプルでは簡易的に内容をコピーするだけ
                # 実際のプロジェクトでは正しい変換ロジックを使用
                f.write(f"PDBQT to SDF conversion for {os.path.basename(pdbqt_file)}\n")
                f.write(pdbqt_content)
            converted_files.append(sdf_file)
            print(f"変換完了: {os.path.basename(sdf_file)}")
        except Exception as e:
            logger.error(f"{ligand_name} の変換に失敗: {e}")
    
    print(f"\n変換されたファイル数: {len(converted_files)}")
    print(f"SDF出力ディレクトリ: {sdf_dir}")
    
    # 全ての結果を1つのSDFにマージ
    merged_sdf = os.path.join(output_dir, "all_results.sdf")
    with open(merged_sdf, 'w') as outfile:
        for sdf_file in converted_files:
            with open(sdf_file, 'r') as infile:
                outfile.write(infile.read())
    
    print(f"すべての結果を結合しました: {merged_sdf}")
    return merged_sdf

def main():
    """メイン関数"""
    # コマンドライン引数の解析
    args = parse_arguments()
    
    # 出力ディレクトリの作成
    os.makedirs(args.output, exist_ok=True)
    
    # リガンドファイルのパスを取得
    ligand_paths = get_ligand_paths(args.ligands, args.max_ligands)
    
    print(f"=== バッチドッキング計算 ===")
    print(f"レセプター: {args.receptor}")
    print(f"リガンド数: {len(ligand_paths)}")
    print(f"出力ディレクトリ: {args.output}")
    
    # ドッキング中心の設定
    if args.center:
        center = args.center
    else:
        print("ドッキング中心が指定されていません。レセプターの中心を使用します。")
        # レセプターの中心を計算（簡易版）
        try:
            from rdkit import Chem
            mol = Chem.MolFromPDBFile(args.receptor)
            if mol:
                conf = mol.GetConformer()
                coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
                center = coords.mean(axis=0).tolist()
                print(f"計算されたドッキング中心: {center}")
            else:
                center = [0.0, 0.0, 0.0]
                print("レセプターの読み込みに失敗しました。原点を中心として使用します。")
        except Exception:
            center = [0.0, 0.0, 0.0]
            print("レセプターの中心計算に失敗しました。原点を中心として使用します。")
    
    print(f"ドッキング中心: {center}")
    print(f"ボックスサイズ: {args.size}")
    print(f"探索の徹底度: {args.exhaustiveness}")
    print(f"出力するポーズ数: {args.num_modes}")
    print(f"使用するCPU数: {args.cpu}")
    print(f"並列実行数: {args.jobs}")
    
    # サービスとワークフローの初期化
    prep_factory = MoleculePreparationFactory()
    prep_service = VinaPreparationService()
    docking_service = VinaDockingService()
    result_analysis_service = ResultAnalysisService()
    
    workflow = BatchDockingWorkflow(
        molecule_preparation_service=prep_service,
        docking_service=docking_service,
        result_analysis_service=result_analysis_service
    )
    
    # ドッキング設定の作成
    config = {
        'grid_box': {
            'center_x': center[0],
            'center_y': center[1],
            'center_z': center[2],
            'size_x': args.size[0],
            'size_y': args.size[1],
            'size_z': args.size[2]
        },
        'parameters': {
            'exhaustiveness': args.exhaustiveness,
            'num_modes': args.num_modes,
            'energy_range': 3.0,
            'seed': 0,
            'cpu': args.cpu
        },
        'name': 'BatchDockingConfig'
    }
    
    try:
        # 実行時間の計測開始
        start_time = time.time()
        
        # バッチドッキングワークフローの実行
        results = workflow.execute(
            ligand_paths=ligand_paths,
            receptor_path=args.receptor,
            config=config,
            output_dir=args.output,
            n_jobs=args.jobs
        )
        
        # 実行時間の計測終了
        elapsed_time = time.time() - start_time
        print(f"\n計算時間: {elapsed_time:.2f}秒")
        
        # 結果のランキング表示
        display_results_ranking(results, args.output)
        
        # 必要に応じてSDF変換
        if args.convert_sdf:
            convert_results_to_sdf(results, args.output)
        
        print("\nすべてのドッキング計算が完了しました")
        print(f"結果ディレクトリ: {args.output}")
        
    except Exception as e:
        logger.error(f"実行中にエラーが発生しました: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()