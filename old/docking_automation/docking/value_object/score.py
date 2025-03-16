from dataclasses import dataclass, field
from enum import Enum, auto
from typing import Optional, Dict, Any, List, Union


class ScoreType(Enum):
    """ドッキングスコアの種類を表す列挙型"""
    BINDING_AFFINITY = auto()  # 結合親和性（kcal/mol）
    BINDING_FREE_ENERGY = auto()  # 結合自由エネルギー（kcal/mol）
    INTERACTION_ENERGY = auto()  # 相互作用エネルギー
    SHAPE_COMPLEMENTARITY = auto()  # 形状相補性
    ELECTROSTATIC = auto()  # 静電相互作用
    HYDROGEN_BOND = auto()  # 水素結合
    HYDROPHOBIC = auto()  # 疎水性相互作用
    VDW = auto()  # ファンデルワールス相互作用
    DESOLVATION = auto()  # 脱溶媒化
    ENTROPY = auto()  # エントロピー
    CUSTOM = auto()  # カスタムスコア
    
    @classmethod
    def from_string(cls, name: str) -> 'ScoreType':
        """文字列からスコアタイプを取得"""
        name_map = {
            'affinity': cls.BINDING_AFFINITY,
            'binding_affinity': cls.BINDING_AFFINITY,
            'energy': cls.BINDING_FREE_ENERGY,
            'binding_free_energy': cls.BINDING_FREE_ENERGY,
            'interaction': cls.INTERACTION_ENERGY,
            'interaction_energy': cls.INTERACTION_ENERGY,
            'shape': cls.SHAPE_COMPLEMENTARITY,
            'shape_complementarity': cls.SHAPE_COMPLEMENTARITY,
            'electrostatic': cls.ELECTROSTATIC,
            'hydrogen_bond': cls.HYDROGEN_BOND,
            'hydrophobic': cls.HYDROPHOBIC,
            'vdw': cls.VDW,
            'desolvation': cls.DESOLVATION,
            'entropy': cls.ENTROPY,
        }
        return name_map.get(name.lower(), cls.CUSTOM)


@dataclass(frozen=True)
class Score:
    """ドッキングスコアを表す値オブジェクト"""
    value: float
    score_type: ScoreType
    name: str
    unit: Optional[str] = None
    uncertainty: Optional[float] = None
    better_is_lower: bool = True  # スコアが低いほど良いかどうか
    metadata: Dict[str, Any] = field(default_factory=dict)
    
    def __post_init__(self) -> None:
        """初期化後の処理"""
        # 追加の検証があれば実装
    
    def is_better_than(self, other: 'Score') -> bool:
        """このスコアが他のスコアより良いかどうかを判定
        
        同じタイプのスコアでない場合は比較できないため、Falseを返す
        """
        if self.score_type != other.score_type:
            return False
        
        if self.better_is_lower:
            return self.value < other.value
        else:
            return self.value > other.value
    
    def get_normalized_value(self, min_val: float, max_val: float) -> float:
        """スコア値を正規化（0～1の範囲に変換）
        
        Args:
            min_val: 最小値
            max_val: 最大値
            
        Returns:
            正規化されたスコア値（0～1）
            
        Note:
            better_is_lower=Trueの場合、最小値が1、最大値が0に対応
            better_is_lower=Falseの場合、最小値が0、最大値が1に対応
        """
        if min_val >= max_val:
            raise ValueError("min_val must be less than max_val")
        
        # 値を0～1の範囲に正規化
        normalized = (self.value - min_val) / (max_val - min_val)
        
        # スコアが低いほど良い場合は反転
        if self.better_is_lower:
            normalized = 1.0 - normalized
        
        return max(0.0, min(1.0, normalized))
    
    def format_with_unit(self) -> str:
        """単位付きでスコアを文字列表現"""
        if self.unit:
            return f"{self.value:.2f} {self.unit}"
        return f"{self.value:.2f}"
    
    def to_dict(self) -> Dict[str, Any]:
        """辞書形式に変換"""
        result = {
            'value': self.value,
            'score_type': self.score_type.name,
            'name': self.name,
            'better_is_lower': self.better_is_lower,
        }
        
        if self.unit:
            result['unit'] = self.unit
            
        if self.uncertainty is not None:
            result['uncertainty'] = self.uncertainty
            
        if self.metadata:
            result['metadata'] = self.metadata
            
        return result
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> 'Score':
        """辞書からScoreを作成"""
        score_type = ScoreType[data['score_type']] if 'score_type' in data else ScoreType.CUSTOM
        
        return cls(
            value=float(data['value']),
            score_type=score_type,
            name=data['name'],
            unit=data.get('unit'),
            uncertainty=data.get('uncertainty'),
            better_is_lower=data.get('better_is_lower', True),
            metadata=data.get('metadata', {})
        )
    
    @classmethod
    def create_binding_affinity(cls, value: float, uncertainty: Optional[float] = None) -> 'Score':
        """結合親和性スコアを作成
        
        Args:
            value: 結合親和性（kcal/mol）
            uncertainty: 不確かさ
            
        Returns:
            結合親和性を表すScoreオブジェクト
        """
        return cls(
            value=value,
            score_type=ScoreType.BINDING_AFFINITY,
            name="Binding Affinity",
            unit="kcal/mol",
            uncertainty=uncertainty,
            better_is_lower=True
        )


@dataclass(frozen=True)
class ScoreSet:
    """複数のドッキングスコアをまとめた値オブジェクト"""
    scores: List[Score]
    
    def get_score(self, score_type: ScoreType) -> Optional[Score]:
        """指定したタイプのスコアを取得"""
        for score in self.scores:
            if score.score_type == score_type:
                return score
        return None
    
    def get_score_by_name(self, name: str) -> Optional[Score]:
        """指定した名前のスコアを取得"""
        for score in self.scores:
            if score.name.lower() == name.lower():
                return score
        return None
    
    def get_best_score(self) -> Optional[Score]:
        """主要スコア（通常は結合親和性）を取得"""
        # 結合親和性のスコアを探す
        binding_affinity = self.get_score(ScoreType.BINDING_AFFINITY)
        if binding_affinity:
            return binding_affinity
        
        # 結合自由エネルギーのスコアを探す
        binding_energy = self.get_score(ScoreType.BINDING_FREE_ENERGY)
        if binding_energy:
            return binding_energy
        
        # スコアがない場合はNoneを返す
        if not self.scores:
            return None
        
        # 最初のスコアを返す
        return self.scores[0]
    
    def to_dict(self) -> Dict[str, float]:
        """スコアを名前:値のマップに変換"""
        return {score.name: score.value for score in self.scores}