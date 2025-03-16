import os
import uuid
import time
import tempfile
import shutil
import logging
import random
from typing import Optional, List, Dict, Any

from docking_automation.docking.service.docking_service import DockingService
from docking_automation.docking.entity.docking_task import DockingTask, TaskStatus
from docking_automation.docking.entity.docking_result import DockingResult
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.docking_parameter import DockingParameters, DockingParameter, ParameterType
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.pose import Pose
from docking_automation.docking.value_object.score import Score, ScoreType
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond

# ロガーの設定
logger = logging.getLogger(__name__)


class MockDockingService(DockingService):
    """ドッキング計算サービスのモック実装
    
    外部ツールを使わずに、ダミーデータを返すモック実装。
    """
    
    # デフォルトのドッキング設定パラメータ
    DEFAULT_PARAMS = {
        "exhaustiveness": 8,
        "num_modes": 9,
        "energy_range": 3,
        "seed": 0,
        "cpu": 1
    }
    
    def __init__(self) -> None:
        """コンストラクタ"""
        logger.info("Initializing MockDockingService")
        
        # 一時ディレクトリの作成
        self._temp_dir = tempfile.mkdtemp(prefix="docking_mock_")
        logger.info(f"Created temporary directory: {self._temp_dir}")
        
        # 結果キャッシュ
        self._result_cache: Dict[str, DockingResult] = {}
    
    def execute(self, task: DockingTask) -> DockingResult:
        """ドッキング計算を実行
        
        Args:
            task: ドッキング計算タスク
            
        Returns:
            ドッキング計算結果
            
        Raises:
            ValueError: タスクが無効な場合
            RuntimeError: 計算の実行に失敗した場合
        """
        # タスクの妥当性を検証
        if not self.validate_task(task):
            raise ValueError(f"Invalid docking task: {task.id}")
        
        # キャッシュに結果があれば返す
        cached_result = self.get_result_for_task(task.id)
        if cached_result:
            return cached_result
        
        # タスクをREADYに更新
        if not task.prepare_for_execution():
            raise ValueError(f"Failed to prepare task for execution: {task.id}")
        
        # タスクを実行中にマーク
        task.mark_as_running()
        
        try:
            # 開始時間を記録
            start_time = time.time()
            
            # 処理時間をシミュレート（1〜3秒）
            simulation_time = random.uniform(1.0, 3.0)
            logger.info(f"Simulating docking calculation for {simulation_time:.2f} seconds")
            time.sleep(simulation_time)
            
            # 出力ディレクトリの作成
            output_dir = os.path.join(self._temp_dir, f"task_{task.id}")
            os.makedirs(output_dir, exist_ok=True)
            
            # ダミーの出力ファイルを作成
            output_path = os.path.join(output_dir, "docking_result.pdbqt")
            with open(output_path, "w") as f:
                f.write("REMARK VINA RESULT: 1 -8.7 0.0 0.0\n")
                f.write("REMARK INTER + INTRA:         -10.271\n")
                f.write("REMARK INTER:                 -10.214\n")
                f.write("REMARK INTRA:                  -0.057\n")
                f.write("REMARK CONF INDEPENDENT:       -3.777\n")
                f.write("REMARK CONF DEPENDENT:         -6.494\n")
                f.write("REMARK  Name = MOL1\n")
                f.write("REMARK  0 active torsions:\n")
                f.write("REMARK  status: ('Active', 'In_Complex')\n")
                f.write("ROOT\n")
                f.write("ATOM      1  C1  MOL     1      -1.234   0.975   0.482  1.00  0.00     0.143 C\n")
                f.write("ATOM      2  C2  MOL     1      -2.649   0.544   0.023  1.00  0.00     0.063 C\n")
                f.write("ATOM      3  O1  MOL     1      -3.250   1.530  -0.818  1.00  0.00    -0.393 OA\n")
                f.write("ENDROOT\n")
                f.write("TORSDOF 0\n")
            
            log_path = os.path.join(output_dir, "docking_log.txt")
            with open(log_path, "w") as f:
                f.write("#################################################################\n")
                f.write("# Mock Docking Log File\n")
                f.write("#################################################################\n")
                f.write("\n")
                f.write("Scoring function : vinardo\n")
                f.write("Rigid receptor: /path/to/receptor.pdbqt\n")
                f.write("Ligand: /path/to/ligand.pdbqt\n")
                f.write("Center: X=0.0 Y=0.0 Z=0.0\n")
                f.write("Size: X=20.0 Y=20.0 Z=20.0\n")
                f.write("Exhaustiveness: 8\n")
                f.write("CPU: 1\n")
                f.write("\n")
                f.write("mode |   affinity | dist from best mode\n")
                f.write("     | (kcal/mol) | rmsd l.b.| rmsd u.b.\n")
                f.write("-----+------------+----------+----------\n")
                f.write("   1       -8.7      0.000      0.000\n")
                f.write("   2       -8.4      1.732      2.525\n")
                f.write("   3       -7.9      1.846      3.125\n")
            
            # 終了時間を記録
            end_time = time.time()
            execution_time = end_time - start_time
            
            # ポーズを作成
            poses = self._create_mock_poses(3)  # 3つのポーズを生成
            
            # スコアを作成
            scores = [
                Score.create_binding_affinity(-8.7),  # 1位のスコア
                Score.create_binding_affinity(-8.4),  # 2位のスコア
                Score.create_binding_affinity(-7.9)   # 3位のスコア
            ]
            
            # 使用ツール情報
            tool_info = {
                "tool": "Mock Docking Tool",
                "version": "1.0.0",
                "command": "mock_docking --receptor receptor.pdbqt --ligand ligand.pdbqt --center 0,0,0 --size 20,20,20"
            }
            
            # 結果オブジェクトの作成
            result = DockingResult(
                id=str(uuid.uuid4()),
                task=task,
                poses=poses,
                scores=scores,
                created_at=end_time,
                tool_info=tool_info,
                execution_time=execution_time,
                metadata={
                    "output_path": output_path,
                    "log_path": log_path
                }
            )
            
            # ポーズごとのスコアを設定
            for i, (pose, score) in enumerate(zip(poses, scores)):
                result.add_pose_score(pose.rank, score)
            
            # タスクを完了にマーク
            task.mark_as_completed()
            
            # 結果をキャッシュに保存
            self._result_cache[task.id] = result
            
            return result
            
        except Exception as e:
            # エラー発生時の処理
            task.mark_as_failed(str(e))
            raise RuntimeError(f"Error executing docking task: {str(e)}")
    
    def create_task(
        self, 
        ligand: Ligand, 
        receptor: Receptor, 
        configuration: DockingConfiguration,
        metadata: Optional[Dict[str, Any]] = None
    ) -> DockingTask:
        """ドッキングタスクを作成
        
        Args:
            ligand: ドッキング対象のリガンド
            receptor: ドッキング対象のレセプター
            configuration: ドッキング設定
            metadata: 追加のメタデータ
            
        Returns:
            生成されたドッキングタスク
        """
        # タスクの作成
        task = DockingTask(
            id=str(uuid.uuid4()),
            ligand=ligand,
            receptor=receptor,
            configuration=configuration,
            status=TaskStatus.PENDING,
            metadata=metadata or {}
        )
        
        return task
    
    def validate_task(self, task: DockingTask) -> bool:
        """タスクの妥当性を検証
        
        Args:
            task: 検証するドッキングタスク
            
        Returns:
            タスクが有効かどうか
        """
        # モック実装では常に有効とする
        return True
    
    def get_supported_parameters(self) -> List[str]:
        """このサービスがサポートするパラメータの一覧を取得
        
        Returns:
            サポートされているパラメータ名のリスト
        """
        return list(self.DEFAULT_PARAMS.keys())
    
    def get_default_configuration(self) -> DockingConfiguration:
        """このサービスのデフォルト設定を取得
        
        Returns:
            デフォルトのドッキング設定
        """
        # デフォルトのグリッドボックス
        grid_box = GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
        
        # デフォルトのパラメータ
        parameters = DockingParameters()
        for name, value in self.DEFAULT_PARAMS.items():
            param_type = ParameterType.INTEGER if isinstance(value, int) else ParameterType.FLOAT
            parameters.add(DockingParameter(name=name, value=value, parameter_type=param_type))
        
        return DockingConfiguration(
            grid_box=grid_box,
            parameters=parameters,
            name="Mock Default Configuration"
        )
    
    def cancel_task(self, task: DockingTask) -> bool:
        """実行中のタスクをキャンセル
        
        Args:
            task: キャンセルするタスク
            
        Returns:
            キャンセルが成功したかどうか
        """
        if task.status == TaskStatus.RUNNING:
            task.mark_as_cancelled()
            return True
        return False
    
    def get_result_for_task(self, task_id: str) -> Optional[DockingResult]:
        """タスクIDに対応する結果を取得（キャッシュから）
        
        Args:
            task_id: 取得するタスクのID
            
        Returns:
            タスクの結果（見つからない場合はNone）
        """
        return self._result_cache.get(task_id)
    
    def _create_mock_poses(self, count: int) -> List[Pose]:
        """ダミーのポーズリストを作成
        
        Args:
            count: 作成するポーズの数
            
        Returns:
            ポーズのリスト
        """
        poses = []
        
        for i in range(count):
            # ダミーの原子と結合を作成
            atoms = [
                Atom(atom_id=1, element="C", x=0.0 + i*0.5, y=0.0, z=0.0),
                Atom(atom_id=2, element="O", x=1.0 + i*0.5, y=0.0, z=0.0),
                Atom(atom_id=3, element="N", x=0.0 + i*0.5, y=1.0, z=0.0)
            ]
            
            bonds = [
                Bond(atom1_id=1, atom2_id=2, bond_type="single"),
                Bond(atom1_id=1, atom2_id=3, bond_type="single")
            ]
            
            # 分子構造を作成
            structure = MoleculeStructure(atoms=atoms, bonds=bonds)
            
            # RMSDを計算（1位のポーズからの距離）
            rmsd = 0.0 if i == 0 else random.uniform(1.0, 3.0)
            
            # ポーズを作成
            pose = Pose(
                rank=i + 1,  # 1から始まるランク
                structure=structure,
                rmsd_to_best=rmsd
            )
            
            poses.append(pose)
        
        return poses
    
    def __del__(self) -> None:
        """デストラクタ"""
        # 一時ディレクトリの削除
        try:
            if hasattr(self, '_temp_dir') and os.path.exists(self._temp_dir):
                shutil.rmtree(self._temp_dir)
                logger.info(f"Removed temporary directory: {self._temp_dir}")
        except Exception as e:
            logger.error(f"Error removing temporary directory: {e}")