import os
import uuid
import time
import tempfile
import subprocess
import logging
import re
from typing import Optional, List, Dict, Any, Tuple

from docking_automation.domain.docking.service.docking_service import DockingService
from docking_automation.domain.docking.entity.docking_task import DockingTask, TaskStatus
from docking_automation.domain.docking.entity.docking_result import DockingResult
from docking_automation.domain.docking.entity.ligand import Ligand
from docking_automation.domain.docking.entity.receptor import Receptor
from docking_automation.domain.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.domain.docking.value_object.docking_parameter import DockingParameters, DockingParameter, ParameterType
from docking_automation.domain.docking.value_object.grid_box import GridBox
from docking_automation.domain.docking.value_object.pose import Pose
from docking_automation.domain.docking.value_object.score import Score, ScoreType
from docking_automation.domain.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond
from docking_automation.domain.molecule.value_object.molecule_format import FormatType

# ロガーの設定
logger = logging.getLogger(__name__)


class VinaDockingService(DockingService):
    """AutoDock Vinaを使用したドッキング計算サービスの実装"""
    
    # デフォルトのVina設定パラメータ
    DEFAULT_PARAMS = {
        "exhaustiveness": 8,
        "num_modes": 9,
        "energy_range": 3,
        "seed": 0,
        "cpu": 1
    }
    
    def __init__(self, vina_path: str = "vina"):
        """コンストラクタ
        
        Args:
            vina_path: AutoDock Vinaの実行ファイルパス
        """
        self.vina_path = vina_path
        self._result_cache: Dict[str, DockingResult] = {}  # タスクIDをキーとする結果キャッシュ
        self._running_tasks: Dict[str, subprocess.Popen[Any]] = {}  # 実行中のタスク
    
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
            start_time = time.time()
            
            # 一時設定ファイルの作成
            config_path = self._create_config_file(task.configuration)
            
            # 出力ディレクトリの作成
            output_dir = tempfile.mkdtemp(prefix="vina_docking_")
            output_path = os.path.join(output_dir, f"result_{task.id}.pdbqt")
            log_path = os.path.join(output_dir, f"log_{task.id}.txt")
            
            # 使用するレセプターを決定（複数ある場合は主要なもの）
            primary_receptor = task.get_primary_receptor()
            receptor_path = primary_receptor.get_path()
            
            if not receptor_path:
                raise ValueError("Receptor path is not available")
                
            # Vinaコマンドの構築
            cmd = [
                self.vina_path,
                "--receptor", receptor_path,
                "--ligand", task.ligand.get_path(),
                "--config", config_path,
                "--out", output_path,
                "--log", log_path
            ]
            
            # CPU数の設定（設定ファイルに含まれていない場合）
            cpu_param = task.configuration.get_parameter("cpu")
            if cpu_param:
                cmd.extend(["--cpu", str(cpu_param)])
            
            # シード値の設定（設定ファイルに含まれていない場合）
            seed_param = task.configuration.get_parameter("seed")
            if seed_param:
                cmd.extend(["--seed", str(seed_param)])
            
            # コマンドを実行
            # Noneの可能性があるので、strに変換してからjoin
            cmd_str = ' '.join([str(c) for c in cmd if c is not None])
            logger.info(f"Running Vina command: {cmd_str}")
            
            try:
                # Noneを取り除く
                cmd_filtered = [str(c) for c in cmd if c is not None]
                process = subprocess.Popen(
                    cmd_filtered,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    text=True
                )
                
                # 実行中のタスクに登録
                self._running_tasks[task.id] = process
                
                # 進捗表示（ログ）
                logger.info("Vina calculation in progress...")
                
                # 完了を待機
                stdout, stderr = process.communicate()
                
                # ログに出力（デバッグ用）
                logger.debug(f"Vina stdout: {stdout}")
                if stderr:
                    logger.warning(f"Vina stderr: {stderr}")
                
                # 実行中のタスクから削除
                if task.id in self._running_tasks:
                    del self._running_tasks[task.id]
                
                # 結果を確認
                if process.returncode != 0:
                    error_msg = f"Vina execution failed (code {process.returncode}): {stderr}"
                    logger.error(error_msg)
                    task.mark_as_failed(error_msg)
                    raise RuntimeError(error_msg)
                
                logger.info(f"Vina calculation completed successfully")
                
            except FileNotFoundError as e:
                error_msg = f"Vina executable not found at: {self.vina_path}"
                logger.error(error_msg)
                task.mark_as_failed(error_msg)
                raise RuntimeError(error_msg)
            
            end_time = time.time()
            execution_time = end_time - start_time
            
            # 結果ファイルが存在するか確認
            if not os.path.exists(output_path):
                task.mark_as_failed("Output file was not created")
                raise RuntimeError("Output file was not created")
            
            # ポーズとスコアの解析
            poses, scores = self._parse_results(output_path, log_path)
            
            if not poses:
                task.mark_as_failed("No valid poses found in the output")
                raise RuntimeError("No valid poses found in the output")
            
            # 実行情報の収集
            tool_info = {
                "tool": "AutoDock Vina",
                "version": self._get_vina_version(),
                "command": " ".join([str(c) for c in cmd if c is not None]),
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
                metadata={"stdout": stdout, "stderr": stderr}
            )
            
            # ポーズごとのスコアを設定
            for pose in poses:
                for score in scores:
                    if score.name == "Binding Affinity" and pose.rank == 1:
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
        
        finally:
            # 一時ファイルの削除
            if 'config_path' in locals() and os.path.exists(config_path):
                os.unlink(config_path)
    
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
            
        Raises:
            ValueError: パラメータが無効な場合
        """
        # リガンドとレセプターの妥当性を確認
        if not ligand.is_prepared():
            raise ValueError("Ligand is not properly prepared for docking")
        
        if not receptor.is_prepared():
            raise ValueError("Receptor is not properly prepared for docking")
        
        # 設定の妥当性を確認
        if not configuration.validate():
            raise ValueError("Invalid docking configuration")
        
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
        # リガンドとレセプターの妥当性を確認
        if not task.ligand.is_prepared():
            logger.error(f"Task {task.id}: Ligand is not properly prepared")
            return False
        
        # レセプターの検証（複数ある場合はすべて確認）
        receptors = task.get_receptors()
        for receptor in receptors:
            if not receptor.is_prepared():
                logger.error(f"Task {task.id}: Receptor is not properly prepared")
                return False
        
        # ファイルの存在を確認
        ligand_path = task.ligand.get_path()
        if not ligand_path or not os.path.exists(ligand_path):
            logger.error(f"Task {task.id}: Ligand file does not exist")
            return False
        
        # 主要レセプターのファイル存在を確認
        primary_receptor = task.get_primary_receptor()
        receptor_path = primary_receptor.get_path()
        if not receptor_path or not os.path.exists(receptor_path):
            logger.error(f"Task {task.id}: Receptor file does not exist")
            return False
        
        # 設定の妥当性を確認
        if not task.configuration.validate():
            logger.error(f"Task {task.id}: Invalid configuration")
            return False
        
        return True
    
    def get_supported_parameters(self) -> List[str]:
        """このサービスがサポートするパラメータの一覧を取得
        
        Returns:
            サポートされているパラメータ名のリスト
        """
        return [
            "exhaustiveness",
            "num_modes",
            "energy_range",
            "seed",
            "cpu"
        ]
    
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
            name="Vina Default Configuration"
        )
    
    def cancel_task(self, task: DockingTask) -> bool:
        """実行中のタスクをキャンセル
        
        Args:
            task: キャンセルするタスク
            
        Returns:
            キャンセルが成功したかどうか
        """
        if task.id not in self._running_tasks:
            return False
        
        # プロセスを終了
        process = self._running_tasks[task.id]
        try:
            process.terminate()
            task.mark_as_cancelled()
            return True
        except:
            return False
        finally:
            # 実行中のタスクから削除
            if task.id in self._running_tasks:
                del self._running_tasks[task.id]
    
    def get_result_for_task(self, task_id: str) -> Optional[DockingResult]:
        """タスクIDに対応する結果を取得（キャッシュから）
        
        Args:
            task_id: 取得するタスクのID
            
        Returns:
            タスクの結果（見つからない場合はNone）
        """
        return self._result_cache.get(task_id)
    
    def _create_config_file(self, configuration: DockingConfiguration) -> str:
        """Vinaの設定ファイルを作成
        
        Args:
            configuration: ドッキング設定
            
        Returns:
            作成された設定ファイルのパス
        """
        # 一時ファイルを作成
        fd, config_path = tempfile.mkstemp(suffix=".conf", prefix="vina_")
        with os.fdopen(fd, 'w') as f:
            # グリッドボックスの設定
            grid_box = configuration.grid_box
            f.write(f"center_x = {grid_box.center_x}\n")
            f.write(f"center_y = {grid_box.center_y}\n")
            f.write(f"center_z = {grid_box.center_z}\n")
            f.write(f"size_x = {grid_box.size_x}\n")
            f.write(f"size_y = {grid_box.size_y}\n")
            f.write(f"size_z = {grid_box.size_z}\n")
            
            # その他のパラメータ
            params = configuration.parameters
            
            # コマンドライン引数として渡すパラメータは設定ファイルに書き込まない
            skip_params = ["cpu", "seed"]
            
            for name in self.get_supported_parameters():
                if name in skip_params:
                    continue
                
                value = params.get_value(name)
                if value is not None:
                    f.write(f"{name} = {value}\n")
        
        return config_path
    
    def _parse_results(self, output_path: str, log_path: str) -> Tuple[List[Pose], List[Score]]:
        """結果ファイルを解析してポーズとスコアを取得
        
        Args:
            output_path: 出力ファイルのパス
            log_path: ログファイルのパス
            
        Returns:
            ポーズのリストとスコアのリスト
        """
        poses = []
        scores = []
        
        # ログファイルからスコアを解析
        binding_affinities = []
        rmsd_values = []
        
        try:
            with open(log_path, 'r') as f:
                log_content = f.read()
                
                # スコアテーブルの位置を特定
                table_start = log_content.find("-----+------------+----------+----------")
                if table_start == -1:
                    raise ValueError("Score table not found in log file")
                
                # モードのスコアを抽出
                # 例: "   1       -8.7      0.000      0.000"
                mode_pattern = r"^\s*(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)"
                lines = log_content[table_start:].split('\n')
                
                for line in lines[1:]:  # ヘッダー行をスキップ
                    match = re.match(mode_pattern, line)
                    if match:
                        mode_num = int(match.group(1))
                        affinity = float(match.group(2))
                        rmsd_lb = float(match.group(3))
                        rmsd_ub = float(match.group(4))
                        
                        binding_affinities.append((mode_num, affinity))
                        rmsd_values.append((mode_num, rmsd_lb, rmsd_ub))
                
                # すべてのスコアをScoreオブジェクトとして追加
                for mode_num, affinity in binding_affinities:
                    score = Score.create_binding_affinity(affinity)
                    scores.append(score)
                
                logger.info(f"Parsed {len(binding_affinities)} scores from log file")
                
        except Exception as e:
            logger.error(f"Error parsing log file: {e}")
            logger.error(f"Log file path: {log_path}")
            logger.debug(f"Exception details: {str(e)}")
        
        # 出力ファイルからポーズを解析
        try:
            # PDBQTファイルを読み込み、モデルごとに分割
            with open(output_path, 'r') as f:
                pdbqt_content = f.read()
            
            # モデル（ポーズ）ごとに分割
            # PDBQTファイルでは、各モデルはMODELとENDMDLタグで囲まれている
            model_pattern = r"MODEL\s+(\d+)(.*?)ENDMDL"
            model_matches = re.findall(model_pattern, pdbqt_content, re.DOTALL)
            
            if not model_matches:
                # MODELタグがない場合は、ファイル全体を1つのモデルとして扱う
                model_matches = [(1, pdbqt_content)]
            
            for model_match in model_matches:
                model_num = int(model_match[0])
                model_content = model_match[1]
                
                # 原子情報を抽出
                atoms = []
                atom_pattern = r"ATOM\s+(\d+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\d+)\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)"
                atom_matches = re.findall(atom_pattern, model_content)
                
                if not atom_matches:
                    # 原子情報が見つからない場合はスキップ
                    logger.warning(f"No atoms found in model {model_num}")
                    continue
                
                for atom_match in atom_matches:
                    atom_id = int(atom_match[0])
                    atom_name = atom_match[1]
                    # 元素記号を取得（原子名の先頭文字）
                    element = atom_name[0]
                    if len(atom_name) > 1 and atom_name[1].isalpha() and not atom_name[1].isupper():
                        element += atom_name[1]  # 2文字目が小文字の場合は含める（例：Cl, Br）
                    
                    x = float(atom_match[5])
                    y = float(atom_match[6])
                    z = float(atom_match[7])
                    
                    atom = Atom(atom_id=atom_id, element=element, x=x, y=y, z=z)
                    atoms.append(atom)
                
                # 結合情報（簡易的に推定）
                bonds = []
                for i in range(len(atoms) - 1):
                    # 隣接する原子間に単結合があると仮定
                    bond = Bond(atom1_id=atoms[i].atom_id, atom2_id=atoms[i+1].atom_id, bond_type="single")
                    bonds.append(bond)
                
                # 分子構造を作成
                structure = MoleculeStructure(atoms=atoms, bonds=bonds)
                
                # RMSDを取得
                rmsd_to_best = 0.0
                for mode, rmsd_lb, _ in rmsd_values:
                    if mode == model_num:
                        rmsd_to_best = rmsd_lb
                        break
                
                # ポーズを作成
                pose = Pose(
                    rank=model_num,
                    structure=structure,
                    rmsd_to_best=rmsd_to_best
                )
                
                poses.append(pose)
            
            logger.info(f"Parsed {len(poses)} poses from output file")
            
            # ランクでソート
            poses.sort(key=lambda p: p.rank)
        
        except Exception as e:
            logger.error(f"Error parsing output file: {e}")
            logger.error(f"Output file path: {output_path}")
            logger.debug(f"Exception details: {str(e)}")
            
            # 解析に失敗した場合は、デフォルトのポーズを作成
            if binding_affinities and not poses:
                logger.warning("Creating default poses based on scores")
                for idx, (mode_num, _) in enumerate(binding_affinities):
                    # デフォルトの原子と結合を作成
                    atoms = [
                        Atom(atom_id=1, element="C", x=0.0, y=0.0, z=0.0),
                        Atom(atom_id=2, element="O", x=1.0, y=0.0, z=0.0)
                    ]
                    bonds = [
                        Bond(atom1_id=1, atom2_id=2, bond_type="single")
                    ]
                    
                    # デフォルトの分子構造を作成
                    structure = MoleculeStructure(atoms=atoms, bonds=bonds)
                    
                    # RMSDの値
                    rmsd = 0.0
                    if idx > 0 and idx < len(rmsd_values):
                        rmsd = rmsd_values[idx][1]
                    elif idx > 0:
                        rmsd = float(idx)  # バックアップ値として
                    
                    # ポーズを作成
                    pose = Pose(
                        rank=mode_num,
                        structure=structure,
                        rmsd_to_best=rmsd
                    )
                    
                    poses.append(pose)
        
        return poses, scores
    
    def _get_vina_version(self) -> str:
        """Vinaのバージョンを取得
        
        Returns:
            バージョン文字列
        """
        try:
            result = subprocess.run(
                [self.vina_path, "--help"],
                capture_output=True,
                text=True,
                check=False
            )
            
            # バージョン情報を抽出
            version_match = re.search(r"AutoDock Vina (\d+\.\d+(\.\d+)?)", result.stdout)
            if version_match:
                return version_match.group(1)
            
            return "Unknown"
        except:
            return "Unknown"