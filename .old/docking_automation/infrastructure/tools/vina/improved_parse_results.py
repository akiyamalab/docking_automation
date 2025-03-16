import re
import logging
from typing import List, Tuple, Optional, Any

from docking_automation.docking.value_object.pose import Pose
from docking_automation.docking.value_object.score import Score
from docking_automation.molecule.value_object.molecule_structure import MoleculeStructure, Atom, Bond

# ロガーの設定
logger = logging.getLogger(__name__)


def _parse_results(self: Any, output_path: str, log_path: str) -> Tuple[List[Pose], List[Score]]:
    """結果ファイルを解析してポーズとスコアを取得
    
    Args:
        output_path: 出力ファイルのパス
        log_path: ログファイルのパス
        
    Returns:
        ポーズのリストとスコアのリスト
    """
    poses: List[Pose] = []
    scores: List[Score] = []
    
    # ログファイルからスコアを解析
    binding_affinities: List[Tuple[int, float]] = []
    rmsd_values: List[Tuple[int, float, float]] = []
    
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
            atoms: List[Atom] = []
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
            bonds: List[Bond] = []
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