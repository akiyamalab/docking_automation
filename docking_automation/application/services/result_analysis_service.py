"""ドッキング結果解析サービス

このモジュールは、ドッキング結果の解析と可視化のためのサービスを提供します。
"""

import os
import logging
import statistics
from typing import List, Dict, Any, Optional, Union
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

from docking_automation.docking.entity.docking_result import DockingResult


class ResultAnalysisService:
    """ドッキング結果解析サービス
    
    ドッキング結果の解析と可視化を行うためのサービスです。
    主な機能：
    - スコアの統計解析
    - スコア分布のヒストグラム生成
    - 相互作用の解析
    """
    
    def __init__(self, logger: Optional[logging.Logger] = None):
        """コンストラクタ
        
        Args:
            logger: ロガー（指定しない場合はモジュールのロガーを使用）
        """
        self.logger = logger or logging.getLogger(__name__)
    
    def generate_score_histogram(
        self, 
        results: List[DockingResult],
        output_path: Union[str, Path],
        bins: int = 20,
        title: str = 'ドッキングスコア分布',
        figsize: tuple = (10, 6),
        color: str = 'skyblue',
        dpi: int = 300
    ) -> str:
        """ドッキングスコアのヒストグラムを生成する
        
        Args:
            results: ドッキング結果のリスト
            output_path: 出力ファイルのパス
            bins: ヒストグラムのビン数
            title: グラフのタイトル
            figsize: 図のサイズ (width, height)
            color: ヒストグラムの色
            dpi: 解像度 (dots per inch)
            
        Returns:
            出力ファイルのパス
            
        Raises:
            ValueError: 結果が空の場合
            IOError: ファイル出力に失敗した場合
        """
        if not results:
            raise ValueError("結果リストが空です")
        
        # 各結果からベストスコアを抽出
        scores = []
        for result in results:
            best_score = result.get_best_score()
            if best_score:
                scores.append(best_score.value)
        
        if not scores:
            raise ValueError("有効なスコアが見つかりません")
        
        # 統計情報の計算
        stats = self.calculate_score_statistics(scores)
        
        # ヒストグラムの生成
        plt.figure(figsize=figsize)
        plt.hist(scores, bins=bins, alpha=0.7, color=color)
        plt.xlabel('ドッキングスコア (kcal/mol)', fontsize=12)
        plt.ylabel('化合物数', fontsize=12)
        plt.title(title, fontsize=14)
        plt.grid(True, alpha=0.3)
        
        # 統計情報をグラフに追加
        info_text = (
            f"サンプル数: {stats['count']}\n"
            f"平均値: {stats['mean']:.2f}\n"
            f"最小値: {stats['min']:.2f}\n"
            f"最大値: {stats['max']:.2f}\n"
            f"標準偏差: {stats['std']:.2f}"
        )
        plt.annotate(
            info_text, 
            xy=(0.02, 0.95), 
            xycoords='axes fraction', 
            bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.8)
        )
        
        # 平均値と最良値を縦線で表示
        plt.axvline(x=stats['mean'], color='red', linestyle='--', alpha=0.7, label=f"平均: {stats['mean']:.2f}")
        plt.axvline(x=stats['min'], color='green', linestyle='--', alpha=0.7, label=f"最良: {stats['min']:.2f}")
        plt.legend()
        
        # 保存
        os.makedirs(os.path.dirname(str(output_path)), exist_ok=True)
        plt.savefig(output_path, dpi=dpi, bbox_inches='tight')
        plt.close()
        
        self.logger.info(f"ヒストグラムを保存しました: {output_path}")
        return str(output_path)
    
    def calculate_score_statistics(self, scores: List[float]) -> Dict[str, float]:
        """スコアの統計情報を計算する
        
        Args:
            scores: スコアのリスト
            
        Returns:
            統計情報を含む辞書
        """
        if not scores:
            return {
                'count': 0,
                'mean': float('nan'),
                'min': float('nan'),
                'max': float('nan'),
                'std': float('nan'),
                'median': float('nan')
            }
        
        return {
            'count': len(scores),
            'mean': statistics.mean(scores),
            'min': min(scores),
            'max': max(scores),
            'std': statistics.stdev(scores) if len(scores) > 1 else 0.0,
            'median': statistics.median(scores)
        }
    
    def rank_results_by_score(
        self, 
        results: List[DockingResult],
        reverse: bool = False
    ) -> List[DockingResult]:
        """結果をスコアでランク付けする
        
        Args:
            results: ドッキング結果のリスト
            reverse: True の場合は降順（スコアが高い順）、False の場合は昇順（スコアが低い順）
            
        Returns:
            ランク付けされた結果のリスト
        """
        def get_score(result: DockingResult) -> float:
            best_score = result.get_best_score()
            # スコアがない場合は無限大（またはマイナス無限大）を返して最後にソートされるようにする
            if not best_score:
                return float('inf') if not reverse else float('-inf')
            return best_score.value
        
        return sorted(results, key=get_score, reverse=reverse)
    
    def compare_results(
        self, 
        results_a: List[DockingResult], 
        results_b: List[DockingResult],
        output_path: Optional[Union[str, Path]] = None
    ) -> Dict[str, Any]:
        """2つの結果セットを比較する
        
        Args:
            results_a: 1つ目の結果セット
            results_b: 2つ目の結果セット
            output_path: 出力ファイルのパス（指定した場合は比較グラフを保存）
            
        Returns:
            比較結果の統計情報
        """
        # スコアを抽出（より安全な方法）
        scores_a = []
        for r in results_a:
            score = r.get_best_score()
            if score is not None:
                scores_a.append(score.value)
                
        scores_b = []
        for r in results_b:
            score = r.get_best_score()
            if score is not None:
                scores_b.append(score.value)
        
        # 統計情報を計算
        stats_a = self.calculate_score_statistics(scores_a)
        stats_b = self.calculate_score_statistics(scores_b)
        
        # グラフ生成（オプション）
        if output_path and scores_a and scores_b:
            plt.figure(figsize=(12, 8))
            
            # ヒストグラム
            plt.subplot(2, 1, 1)
            plt.hist(scores_a, bins=20, alpha=0.5, label='セットA', color='blue')
            plt.hist(scores_b, bins=20, alpha=0.5, label='セットB', color='orange')
            plt.xlabel('ドッキングスコア (kcal/mol)')
            plt.ylabel('頻度')
            plt.title('スコア分布の比較')
            plt.legend()
            plt.grid(True, alpha=0.3)
            
            # ボックスプロット
            plt.subplot(2, 1, 2)
            plt.boxplot([scores_a, scores_b], labels=['セットA', 'セットB'])
            plt.ylabel('ドッキングスコア (kcal/mol)')
            plt.title('統計分布の比較')
            plt.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
        
        # 比較結果をまとめる
        return {
            'set_a': stats_a,
            'set_b': stats_b,
            'mean_diff': stats_a['mean'] - stats_b['mean'],
            'min_diff': stats_a['min'] - stats_b['min'],
            'count_a': len(scores_a),
            'count_b': len(scores_b),
        }