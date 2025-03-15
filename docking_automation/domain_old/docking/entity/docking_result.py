import uuid
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Tuple
import time

from docking_automation.domain.docking.entity.docking_task import DockingTask
from docking_automation.domain.docking.value_object.pose import Pose, PoseEnsemble
from docking_automation.domain.docking.value_object.score import Score, ScoreSet
from docking_automation.domain.docking.value_object.docking_parameter import DockingParameters


@dataclass
class DockingResult:
    """ドッキング計算の結果を表すエンティティ"""
    
    id: str
    task: DockingTask
    poses: List[Pose]
    scores: List[Score]
    created_at: float  # Unix timestamp
    tool_info: Dict[str, str] = field(default_factory=dict)  # 使用したツールの情報
    execution_time: Optional[float] = None  # 実行時間（秒）
    metadata: Dict[str, Any] = field(default_factory=dict)
    pose_scores: Dict[int, List[Score]] = field(default_factory=dict)  # ポーズごとのスコア
    
    def __post_init__(self) -> None:
        """初期化後の処理"""
        # IDが指定されていない場合はUUIDを生成
        if self.id == "":
            object.__setattr__(self, "id", uuid.uuid4().hex)
        
        # 作成日時が指定されていない場合は現在時刻を設定
        if self.created_at == 0.0:
            object.__setattr__(self, "created_at", time.time())
    
    def get_best_pose(self) -> Optional[Pose]:
        """最良のポーズを取得（最もスコアが良いポーズ）"""
        if not self.poses:
            return None
        
        # ランクでソートして最初のポーズを返す
        return sorted(self.poses, key=lambda p: p.rank)[0]
    
    def get_pose_ensemble(self) -> PoseEnsemble:
        """すべてのポーズをPoseEnsembleとして取得"""
        return PoseEnsemble(poses=self.poses)
    
    def get_score_set(self) -> ScoreSet:
        """すべてのスコアをScoreSetとして取得"""
        return ScoreSet(scores=self.scores)
    
    def get_best_score(self) -> Optional[Score]:
        """最良のスコアを取得（通常は結合親和性）"""
        return self.get_score_set().get_best_score()
    
    def get_scores_for_pose(self, pose_rank: int) -> List[Score]:
        """指定したランクのポーズに関連するスコアを取得"""
        return self.pose_scores.get(pose_rank, [])
    
    def add_pose_score(self, pose_rank: int, score: Score) -> None:
        """ポーズに関連するスコアを追加"""
        if pose_rank not in self.pose_scores:
            self.pose_scores[pose_rank] = []
        
        self.pose_scores[pose_rank].append(score)
    
    def get_pose_by_rank(self, rank: int) -> Optional[Pose]:
        """指定されたランクのポーズを取得"""
        for pose in self.poses:
            if pose.rank == rank:
                return pose
        return None
    
    def get_ranked_poses(self) -> List[Tuple[Pose, List[Score]]]:
        """ポーズとそのスコアのペアをランク順に取得"""
        result = []
        for pose in sorted(self.poses, key=lambda p: p.rank):
            scores = self.get_scores_for_pose(pose.rank)
            result.append((pose, scores))
        return result
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        result = {
            'id': self.id,
            'task_id': self.task.id,
            'created_at': self.created_at,
            'execution_time': self.execution_time,
            'tool_info': self.tool_info,
            'metadata': self.metadata,
            'scores': [score.to_dict() for score in self.scores],
            'best_score': None,
        }
        
        # bestスコアがあれば追加
        best_score = self.get_best_score()
        if best_score:
            result['best_score'] = best_score.to_dict()
        
        # ポーズの情報も含める（必要に応じて）
        if self.poses:
            result['poses'] = {
                'count': len(self.poses),
                'best_pose_rank': None,
            }
            best_pose = self.get_best_pose()
            if best_pose:
                poses_dict = result['poses']
                if isinstance(poses_dict, dict):
                    poses_dict['best_pose_rank'] = best_pose.rank
        
        return result
    
    @classmethod
    def create_empty_result(cls, task: DockingTask) -> 'DockingResult':
        """空の結果オブジェクトを作成（失敗時など）"""
        return cls(
            id=uuid.uuid4().hex,
            task=task,
            poses=[],
            scores=[],
            created_at=time.time(),
            execution_time=0.0,
            tool_info={'status': 'failed'},
            metadata={'error': 'No results generated'}
        )
    
    def __str__(self) -> str:
        """文字列表現"""
        best_score = self.get_best_score()
        score_str = f"{best_score.value:.2f} {best_score.unit}" if best_score else "N/A"
        return f"DockingResult(id={self.id}, poses={len(self.poses)}, best_score={score_str})"
