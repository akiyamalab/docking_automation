"""バッチドッキング計算ワークフロー

このモジュールは、複数の化合物に対するバッチドッキング計算ワークフローを提供します。
"""

import os
import uuid
import logging
from typing import Dict, Any, Optional, List, Union, TYPE_CHECKING, cast
from pathlib import Path
import concurrent.futures

if TYPE_CHECKING:
    import pandas as pd
else:
    import pandas as pd

# 分子モジュール
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.service.molecule_preparation_service import MoleculePreparationService

# ドッキングモジュール
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.entity.docking_task import DockingTask
from docking_automation.docking.entity.docking_result import DockingResult
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.docking_parameter import DockingParameter, ParameterType as DPParameterType
from docking_automation.docking.value_object.docking_parameters import DockingParameters
# あえて別の名前でインポートを避ける
from docking_automation.docking.service.docking_service import DockingService

# 結果解析モジュール（新規追加）
from docking_automation.application.services.result_analysis_service import ResultAnalysisService

# ロガーの設定
logger = logging.getLogger(__name__)


class BatchDockingWorkflow:
    """複数の化合物に対するバッチドッキング計算ワークフロー
    
    複数の化合物ファイルと単一のタンパク質のドッキング計算を行うワークフロー。
    1. 各化合物とタンパク質をロードし、準備する
    2. 並行してドッキング計算を実行する
    3. 結果を統合して返す
    """
    
    def __init__(
        self,
        molecule_preparation_service: MoleculePreparationService,
        docking_service: DockingService,
        result_analysis_service: Optional[ResultAnalysisService] = None
    ):
        """コンストラクタ
        
        Args:
            molecule_preparation_service: 分子準備サービス
            docking_service: ドッキング計算サービス
            result_analysis_service: 結果解析サービス（オプション）
        """
        self.molecule_preparation_service = molecule_preparation_service
        self.docking_service = docking_service
        self.result_analysis_service = result_analysis_service
    
    def execute(
        self,
        ligand_paths: List[Union[str, Path]],
        receptor_path: Union[str, Path],
        config: Optional[Dict[str, Any]] = None,
        output_dir: Optional[Union[str, Path]] = None,
        n_jobs: int = 1
    ) -> List[DockingResult]:
        """バッチドッキング計算を実行する
        
        Args:
            ligand_paths: リガンドファイルのパスのリスト
            receptor_path: レセプターファイルのパス
            config: ドッキング設定（Noneの場合はデフォルト設定を使用）
            output_dir: 出力ディレクトリ（Noneの場合は一時ディレクトリを使用）
            n_jobs: 並列実行数（デフォルト: 1）
            
        Returns:
            ドッキング計算結果のリスト
            
        Raises:
            ValueError: パラメータが無効な場合
            RuntimeError: 処理に失敗した場合
        """
        try:
            # 入力チェック
            if not ligand_paths:
                raise ValueError("リガンドファイルのパスが指定されていません")
            
            for ligand_path in ligand_paths:
                if not os.path.exists(str(ligand_path)):
                    raise ValueError(f"リガンドファイルが見つかりません: {ligand_path}")
            
            if not os.path.exists(str(receptor_path)):
                raise ValueError(f"レセプターファイルが見つかりません: {receptor_path}")
            
            # 出力ディレクトリの作成
            if output_dir and not os.path.exists(str(output_dir)):
                os.makedirs(str(output_dir))
            
            # 1. タンパク質の準備（全リガンドで共通）
            logger.info(f"レセプターの準備: {receptor_path}")
            protein = self._load_protein(str(receptor_path))
            prepared_protein = self.molecule_preparation_service.prepare_receptor(protein)
            receptor = Receptor(protein=prepared_protein)
            
            # 2. ドッキング設定の生成（全リガンドで共通）
            docking_configuration = self._create_configuration(config, receptor)
            
            # 3. 並列ドッキング実行
            results = self._run_parallel_docking(
                ligand_paths, 
                receptor, 
                docking_configuration, 
                output_dir, 
                n_jobs
            )
            
            # 4. バッチ結果のサマリー作成
            if output_dir:
                self._save_batch_summary(results, str(output_dir))
            
            return results
            
        except Exception as e:
            logger.error(f"BatchDockingWorkflowでエラーが発生: {e}")
            raise RuntimeError(f"バッチドッキングワークフローでエラーが発生: {str(e)}")
    
    def _create_configuration(
        self, 
        config: Optional[Dict[str, Any]], 
        receptor: Receptor
    ) -> DockingConfiguration:
        """ドッキング設定を生成する
        
        Args:
            config: ユーザー指定の設定
            receptor: レセプター
            
        Returns:
            ドッキング設定
        """
        if not config:
            # デフォルト設定を使用
            docking_configuration = self.docking_service.get_default_configuration()
            
            # レセプターのグリッドボックスを使用（存在する場合）
            if receptor.has_active_site_defined():
                docking_configuration = docking_configuration.with_grid_box(
                    receptor.get_suggested_grid_box()
                )
            return docking_configuration
        else:
            # 設定から構築
            grid_config = config.get('grid_box', {})
            grid_box = GridBox(
                center_x=float(grid_config.get('center_x', 0.0)),
                center_y=float(grid_config.get('center_y', 0.0)),
                center_z=float(grid_config.get('center_z', 0.0)),
                size_x=float(grid_config.get('size_x', 20.0)),
                size_y=float(grid_config.get('size_y', 20.0)),
                size_z=float(grid_config.get('size_z', 20.0))
            )
            
            # パラメータの設定
            params_config = config.get('parameters', {})
            params = DockingParameters()
            for name, value in params_config.items():
                # mypy向け型ヒント
                if isinstance(value, int):
                    pt = cast(DPParameterType, DPParameterType.INTEGER)
                else:
                    pt = cast(DPParameterType, DPParameterType.FLOAT)
                params.add(DockingParameter(name=name, value=value, parameter_type=pt))
            
            # 設定を作成して返す
            return DockingConfiguration(
                grid_box=grid_box,
                parameters=cast(Any, params),
                name=config.get('name', 'バッチドッキング設定')
            )
    
    def _run_parallel_docking(
        self,
        ligand_paths: List[Union[str, Path]],
        receptor: Receptor,
        docking_configuration: DockingConfiguration,
        output_dir: Optional[Union[str, Path]],
        n_jobs: int
    ) -> List[DockingResult]:
        """複数のリガンドに対して並列ドッキングを実行する
        
        Args:
            ligand_paths: リガンドファイルのパスのリスト
            receptor: 準備済みレセプター
            docking_configuration: ドッキング設定
            output_dir: 出力ディレクトリ
            n_jobs: 並列実行数
            
        Returns:
            ドッキング結果のリスト
        """
        results = []
        
        # 並列実行（ThreadPoolExecutorを使用）
        with concurrent.futures.ThreadPoolExecutor(max_workers=n_jobs) as executor:
            # 各リガンドのドッキングタスクを作成して送信
            future_to_path = {
                executor.submit(
                    self._dock_single_ligand, 
                    ligand_path, 
                    receptor, 
                    docking_configuration,
                    output_dir
                ): ligand_path 
                for ligand_path in ligand_paths
            }
            
            # 完了したタスクから結果を取得
            for future in concurrent.futures.as_completed(future_to_path):
                ligand_path = future_to_path[future]
                try:
                    result = future.result()
                    results.append(result)
                    logger.info(f"リガンド {ligand_path} のドッキングが完了しました")
                except Exception as e:
                    logger.error(f"リガンド {ligand_path} のドッキングでエラーが発生: {e}")
        
        return results
    
    def _dock_single_ligand(
        self,
        ligand_path: Union[str, Path],
        receptor: Receptor,
        docking_configuration: DockingConfiguration,
        output_dir: Optional[Union[str, Path]]
    ) -> DockingResult:
        """単一のリガンドに対してドッキング計算を実行する
        
        Args:
            ligand_path: リガンドファイルのパス
            receptor: 準備済みレセプター
            docking_configuration: ドッキング設定
            output_dir: 出力ディレクトリ
            
        Returns:
            ドッキング結果
        """
        # 1. 化合物の準備
        logger.info(f"リガンドの準備: {ligand_path}")
        compound = self._load_compound(str(ligand_path))
        prepared_compound = self.molecule_preparation_service.prepare_ligand(compound)
        ligand = Ligand(compound=prepared_compound)
        
        # 2. ドッキングタスクの作成
        task = self.docking_service.create_task(
            ligand=ligand,
            receptor=receptor,
            configuration=docking_configuration,
            metadata={
                'ligand_path': str(ligand_path),
                'batch_processing': True
            }
        )
        
        # 3. ドッキング計算の実行
        logger.info(f"ドッキングタスクの実行: {task.id}")
        result = self.docking_service.execute(task)
        
        # 4. 結果の保存（リガンド固有のディレクトリ）
        if output_dir:
            # リガンド名からサブディレクトリ名を生成
            ligand_name = os.path.splitext(os.path.basename(str(ligand_path)))[0]
            ligand_output_dir = os.path.join(str(output_dir), ligand_name)
            os.makedirs(ligand_output_dir, exist_ok=True)
            
            # 結果ファイルの保存
            self._save_result(result, ligand_output_dir)
        
        return result
    
    def _load_compound(self, file_path: str) -> Compound:
        """ファイルから化合物をロードする"""
        # 実装を簡略化（domain削除対応用）
        return Compound(
            id=str(uuid.uuid4()),
            path=file_path
        )
    
    def _load_protein(self, file_path: str) -> Protein:
        """ファイルからタンパク質をロードする"""
        # 実装を簡略化（domain削除対応用）
        return Protein(
            id=str(uuid.uuid4()),
            path=file_path
        )
    
    def _save_result(self, result: DockingResult, output_dir: str) -> None:
        """リガンド単位の結果を保存する"""
        # ベストスコアと実行時間を出力
        best_score = result.get_best_score()
        if best_score:
            score_str = f"{best_score.value:.2f} {best_score.unit}" if best_score.unit else f"{best_score.value:.2f}"
        else:
            score_str = "N/A"
        
        summary_path = os.path.join(output_dir, "docking_summary.txt")
        with open(summary_path, "w") as f:
            f.write(f"Docking Result Summary\n")
            f.write(f"======================\n\n")
            f.write(f"Task ID: {result.task.id}\n")
            f.write(f"Ligand: {result.task.ligand.name}\n")
            f.write(f"Receptor: {result.task.get_primary_receptor().name}\n")
            f.write(f"Poses: {len(result.poses)}\n")
            f.write(f"Best Score: {score_str}\n")
            f.write(f"Execution Time: {result.execution_time:.2f} seconds\n")
    
    def _save_batch_summary(self, results: List[DockingResult], output_dir: str) -> None:
        """バッチドッキング結果のサマリーを保存する
        
        Args:
            results: ドッキング結果のリスト
            output_dir: 出力ディレクトリ
        """
        # CSV形式でランキング結果を保存
        summary_data = []
        
        for result in results:
            ligand_name = os.path.basename(result.task.metadata.get('ligand_path', 'unknown'))
            best_score = result.get_best_score()
            score_value = best_score.value if best_score else float('inf')
            
            summary_data.append({
                'ligand': ligand_name,
                'score': score_value,
                'poses': len(result.poses),
                'execution_time': result.execution_time
            })
        
        # スコア順にソート
        sorted_data = sorted(summary_data, key=lambda x: x['score'])
        
        # DataFrameに変換してCSVとして保存
        df = pd.DataFrame(sorted_data)
        csv_path = os.path.join(output_dir, "ranked_results.csv")
        df.to_csv(csv_path, index=False)
        
        # テキスト形式でも保存
        txt_path = os.path.join(output_dir, "batch_summary.txt")
        with open(txt_path, "w") as f:
            f.write(f"Batch Docking Results\n")
            f.write(f"===================\n\n")
            f.write(f"Total ligands processed: {len(results)}\n\n")
            f.write(f"Ranked Results:\n")
            f.write(f"--------------\n")
            
            for i, data in enumerate(sorted_data):
                f.write(f"{i+1}. {data['ligand']} - Score: {data['score']:.2f}\n")