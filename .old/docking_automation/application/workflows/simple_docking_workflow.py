import os
import uuid
import logging
from typing import Dict, Any, Optional, List, Tuple

# 分子モジュール
from docking_automation.molecule import (
    Compound,
    Protein,
    MoleculePreparationService,
)
from docking_automation.molecule.value_object.molecule_format import FormatType, MoleculeFormat
from docking_automation.molecule.value_object.molecule_structure import Atom, Bond, MoleculeStructure

# ドッキングモジュール
from docking_automation.docking import (
    Ligand,
    Receptor,
    DockingTask,
    DockingResult,
    DockingConfiguration,
    GridBox,
    DockingParameters,
    DockingParameter,
    ParameterType,
    DockingService,
)

# ロガーの設定
logger = logging.getLogger(__name__)


class SimpleDockingWorkflow:
    """シンプルなドッキング計算ワークフロー
    
    単一の化合物とタンパク質のドッキング計算を行うワークフロー。
    1. 化合物とタンパク質をロードし、準備する
    2. ドッキング計算を実行する
    3. 結果を返す
    """
    
    def __init__(
        self,
        molecule_preparation_service: MoleculePreparationService,
        docking_service: DockingService
    ):
        """コンストラクタ
        
        Args:
            molecule_preparation_service: 分子準備サービス
            docking_service: ドッキング計算サービス
        """
        self.molecule_preparation_service = molecule_preparation_service
        self.docking_service = docking_service
    
    def execute(
        self,
        ligand_path: str,
        receptor_path: str,
        config: Optional[Dict[str, Any]] = None,
        output_dir: Optional[str] = None
    ) -> DockingResult:
        """ドッキング計算を実行する
        
        Args:
            ligand_path: リガンドファイルのパス
            receptor_path: レセプターファイルのパス
            config: ドッキング設定（Noneの場合はデフォルト設定を使用）
            output_dir: 出力ディレクトリ（Noneの場合は一時ディレクトリを使用）
            
        Returns:
            ドッキング計算結果
            
        Raises:
            ValueError: パラメータが無効な場合
            RuntimeError: 処理に失敗した場合
        """
        try:
            # 入力ファイルの存在確認
            if not os.path.exists(ligand_path):
                raise ValueError(f"Ligand file not found: {ligand_path}")
            
            if not os.path.exists(receptor_path):
                raise ValueError(f"Receptor file not found: {receptor_path}")
            
            # 出力ディレクトリの作成
            if output_dir and not os.path.exists(output_dir):
                os.makedirs(output_dir)
            
            # 1. 化合物の準備
            logger.info(f"Loading and preparing ligand: {ligand_path}")
            compound = self._load_compound(ligand_path)
            prepared_compound = self.molecule_preparation_service.prepare_ligand(compound)
            
            # 2. タンパク質の準備
            logger.info(f"Loading and preparing receptor: {receptor_path}")
            protein = self._load_protein(receptor_path)
            prepared_protein = self.molecule_preparation_service.prepare_receptor(protein)
            
            # 3. リガンドとレセプターの作成
            ligand = Ligand(compound=prepared_compound)
            receptor = Receptor(protein=prepared_protein)
            
            # 4. ドッキング設定の作成
            # 設定がない場合はレセプターからグリッドボックスを取得
            if not config:
                docking_configuration = self.docking_service.get_default_configuration()
                # レセプターのグリッドボックスを使用（存在する場合）
                if receptor.has_active_site_defined():
                    docking_configuration = docking_configuration.with_grid_box(
                        receptor.get_suggested_grid_box()
                    )
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
                parameters = DockingParameters()
                for name, value in params_config.items():
                    param_type = ParameterType.INTEGER if isinstance(value, int) else ParameterType.FLOAT
                    parameters.add(DockingParameter(name=name, value=value, parameter_type=param_type))
                
                docking_configuration = DockingConfiguration(
                    grid_box=grid_box,
                    parameters=parameters,
                    name=config.get('name', 'Custom Configuration')
                )
            
            # 5. ドッキングタスクの作成
            task = self.docking_service.create_task(
                ligand=ligand,
                receptor=receptor,
                configuration=docking_configuration,
                metadata={
                    'ligand_path': ligand_path,
                    'receptor_path': receptor_path,
                    'output_dir': output_dir
                }
            )
            
            # 6. ドッキング計算の実行
            logger.info(f"Executing docking task: {task.id}")
            result = self.docking_service.execute(task)
            
            # 7. 結果の保存（出力ディレクトリが指定されている場合）
            if output_dir:
                self._save_result(result, output_dir)
            
            return result
            
        except Exception as e:
            logger.error(f"Error in SimpleDockingWorkflow: {e}")
            raise RuntimeError(f"Docking workflow failed: {str(e)}")
    
    def _load_compound(self, file_path: str) -> Compound:
        """ファイルから化合物をロードする"""
        # ファイル形式の推測
        format_type = FormatType.from_path(file_path)
        
        # 実際のアプリケーションでは、ファイルから分子構造を読み込む
        # このサンプル実装では、ダミーの構造を返す
        
        # ダミーの原子と結合を作成
        atoms = [
            Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="O", x=1.0, y=0.0, z=0.0)
        ]
        bonds = [
            Bond(atom1_id=1, atom2_id=2, bond_type="single")
        ]
        
        # 分子構造を作成
        structure = MoleculeStructure(atoms=atoms, bonds=bonds)
        
        # 化合物を作成
        return Compound(
            id=str(uuid.uuid4()),
            structure=structure,
            format=MoleculeFormat(type=format_type),
            path=file_path
        )
    
    def _load_protein(self, file_path: str) -> Protein:
        """ファイルからタンパク質をロードする"""
        # ファイル形式の推測
        format_type = FormatType.from_path(file_path)
        
        # 実際のアプリケーションでは、ファイルから分子構造を読み込む
        # このサンプル実装では、ダミーの構造を返す
        
        # ダミーの原子と結合を作成
        atoms = [
            Atom(atom_id=1, element="N", x=0.0, y=0.0, z=0.0),
            Atom(atom_id=2, element="C", x=1.0, y=0.0, z=0.0),
            Atom(atom_id=3, element="C", x=1.0, y=1.0, z=0.0),
            Atom(atom_id=4, element="O", x=2.0, y=1.0, z=0.0)
        ]
        bonds = [
            Bond(atom1_id=1, atom2_id=2, bond_type="single"),
            Bond(atom1_id=2, atom2_id=3, bond_type="single"),
            Bond(atom1_id=3, atom2_id=4, bond_type="double")
        ]
        
        # 分子構造を作成
        structure = MoleculeStructure(atoms=atoms, bonds=bonds)
        
        # タンパク質を作成
        return Protein(
            id=str(uuid.uuid4()),
            structure=structure,
            format=MoleculeFormat(type=format_type),
            path=file_path,
            chains={"A"}
        )
    
    def _save_result(self, result: DockingResult, output_dir: str) -> None:
        """結果を保存する"""
        # この実装では、単純にベストスコアとポーズ数を出力するだけ
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
            receptor_info = ""
            if result.task.is_multi_receptor():
                receptors = result.task.get_receptors()
                if receptors:
                    primary_name = receptors[0].name if receptors[0].name is not None else "Unnamed"
                    receptor_info = f"{primary_name} (+{len(receptors)-1} more)"
                else:
                    receptor_info = "None"
            else:
                primary_receptor = result.task.get_primary_receptor()
                receptor_info = primary_receptor.name if primary_receptor.name is not None else "Unnamed"
            f.write(f"Receptor: {receptor_info}\n")
            f.write(f"Poses: {len(result.poses)}\n")
            f.write(f"Best Score: {score_str}\n")
            f.write(f"Execution Time: {result.execution_time:.2f} seconds\n")