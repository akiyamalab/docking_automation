from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Tuple

from docking_automation.domain.molecule.value_object.molecule_structure import MoleculeStructure


@dataclass(frozen=True)
class Pose:
    """ドッキング計算の結果として得られる結合ポーズを表す値オブジェクト"""
    
    rank: int  # ポーズのランク（1が最良）
    structure: MoleculeStructure  # ポーズの3D構造
    rmsd_to_best: Optional[float] = None  # 最良ポーズからのRMSD（Å）
    rmsd_to_reference: Optional[float] = None  # リファレンス構造からのRMSD（Å）
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def get_coordinates(self) -> List[Tuple[float, float, float]]:
        """原子の座標リストを取得"""
        return [(atom.x, atom.y, atom.z) for atom in self.structure.atoms]
    
    def get_center(self) -> Tuple[float, float, float]:
        """ポーズの中心座標を取得"""
        return self.structure.get_center_of_mass()
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        result = {
            'rank': self.rank,
            'rmsd_to_best': self.rmsd_to_best,
            'rmsd_to_reference': self.rmsd_to_reference,
            'metadata': self.metadata,
        }
        
        # 構造情報も含める場合はコメントを外す
        # result['structure'] = {
        #     'atoms': [
        #         {
        #             'atom_id': atom.atom_id,
        #             'element': atom.element,
        #             'x': atom.x,
        #             'y': atom.y,
        #             'z': atom.z,
        #         }
        #         for atom in self.structure.atoms
        #     ]
        # }
        
        return result


@dataclass(frozen=True)
class PoseEnsemble:
    """複数の結合ポーズをまとめた値オブジェクト"""
    
    poses: List[Pose]
    
    def get_best_pose(self) -> Optional[Pose]:
        """最良ポーズを取得（ランク1）"""
        if not self.poses:
            return None
        
        # ランクでソートして最初のポーズを返す
        return sorted(self.poses, key=lambda p: p.rank)[0]
    
    def get_pose_by_rank(self, rank: int) -> Optional[Pose]:
        """指定されたランクのポーズを取得"""
        for pose in self.poses:
            if pose.rank == rank:
                return pose
        return None
    
    def get_rmsd_matrix(self) -> List[List[float]]:
        """ポーズ間のRMSDマトリクスを計算
        
        Returns:
            RMSDマトリクス（n x n、nはポーズ数）
        """
        n = len(self.poses)
        matrix = [[0.0 for _ in range(n)] for _ in range(n)]
        
        # 実際のRMSD計算は実装依存のため、ここでは仮の実装
        # （適切な計算メソッドへのプレースホルダー）
        for i in range(n):
            for j in range(i + 1, n):
                rmsd = 0.0
                # 実際にはポーズの座標から計算する
                rmsd_i_val = self.poses[i].rmsd_to_best
                rmsd_j_val = self.poses[j].rmsd_to_best
                
                if rmsd_i_val is not None and rmsd_j_val is not None:
                    # 両方のRMSDが存在する場合のみ計算
                    # この時点でrmsd_i_valとrmsd_j_valはfloatであることが保証されています
                    # 仮の計算：最良ポーズからのRMSDの差の絶対値
                    rmsd = abs(rmsd_i_val - rmsd_j_val)
                
                matrix[i][j] = rmsd
                matrix[j][i] = rmsd
        
        return matrix
    
    def cluster_poses(self, rmsd_threshold: float = 2.0) -> Dict[int, List[Pose]]:
        """ポーズをRMSDに基づいてクラスタリング
        
        Args:
            rmsd_threshold: クラスタリングのRMSDしきい値（Å）
            
        Returns:
            クラスターID:ポーズリストの辞書
        """
        # 実際のクラスタリングアルゴリズムは実装依存
        # ここではプレースホルダーとして、RMSDが最良ポーズからのRMSDに基づく簡易クラスタリングを実装
        
        clusters: Dict[int, List[Pose]] = {}
        cluster_id = 0
        
        for pose in sorted(self.poses, key=lambda p: p.rank):
            # すでに割り当てられたクラスターがあるか確認
            assigned = False
            
            for cid, cluster_poses in clusters.items():
                if not cluster_poses:
                    continue
                    
                representative = cluster_poses[0]
                
                # 両方のRMSDが存在する場合のみ計算
                pose_rmsd_val = pose.rmsd_to_best
                rep_rmsd_val = representative.rmsd_to_best
                
                if pose_rmsd_val is not None and rep_rmsd_val is not None:
                    # この時点でpose_rmsd_valとrep_rmsd_valはfloatであることが保証されています
                    rmsd_diff = abs(pose_rmsd_val - rep_rmsd_val)
                    
                    if rmsd_diff < rmsd_threshold:
                        clusters[cid].append(pose)
                        assigned = True
                        break
            
            # 既存のクラスターに割り当てられなかった場合、新しいクラスターを作成
            if not assigned:
                clusters[cluster_id] = [pose]
                cluster_id += 1
        
        return clusters