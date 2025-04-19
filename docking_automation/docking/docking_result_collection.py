import itertools
from pathlib import Path
from typing import Any, Callable, Dict, Iterator, List, Optional, TypeVar, Union

from sortedcontainers import SortedList  # type: ignore[import-untyped]

from .docking_result import DockingResult

T = TypeVar("T")
DockingResultT = TypeVar("DockingResultT", bound=DockingResult)


class DockingResultCollection:
    """
    複数のDockingResultを管理するコレクションクラス。

    常にドッキングスコアでソートされた状態を維持する。
    このクラスは単純なコレクションとして機能し、DockingResultの集合を管理するためのユーティリティメソッドを提供する。
    大規模なデータセットを扱う場合は、ストリーミング処理やバッチ処理を検討すること。
    """

    # TODO: factoryを使って results list[DockingResult] = [] とする
    def __init__(self, results: Optional[List[DockingResult]] = None, execution_time: Optional[float] = None) -> None:
        """
        DockingResultCollectionオブジェクトを初期化する。

        Args:
            results: 初期結果のリスト（指定しない場合は空のコレクションを作成）
            execution_time: このコレクションが生成される元となった処理の実行時間（秒）（オプション）
        """
        self._results: SortedList[DockingResult] = SortedList(key=lambda result: result.docking_score)
        self.execution_time: Optional[float] = execution_time
        if results:
            for result in results:
                self._results.add(result)

    def add(self, result: DockingResult) -> None:
        """
        結果を追加する。

        Args:
            result: 追加する結果
        """
        self._results.add(result)

    def extend(self, results: Union[List[DockingResult], "DockingResultCollection"]) -> None:
        """
        複数の結果を一度に追加する。

        Args:
            results: 追加する結果のリストまたはDockingResultCollection
        """
        if isinstance(results, DockingResultCollection):
            for result in results:
                self._results.add(result)
        else:
            for result in results:
                self._results.add(result)

    def get_all(self) -> List[DockingResult]:
        """
        すべての結果を取得する。

        Returns:
            結果のリスト
        """
        return list(self._results)

    def get_top(self, n: int) -> List[DockingResult]:
        """
        スコアが上位n件の結果を取得する。

        Args:
            n: 取得する件数

        Returns:
            上位n件の結果のリスト
        """
        # 明示的にリストを作成して型を保証する
        top_n = []
        for i, result in enumerate(self._results):
            if i >= n:
                break
            top_n.append(result)
        return top_n

    def filter(self, condition: Callable[[DockingResult], bool]) -> "DockingResultCollection":
        """
        条件に合致する結果のみを含む新しいコレクションを作成する。

        例えば、特定のタンパク質IDに関連する結果のみを抽出する場合：
        ```python
        protein_results = collection.filter(lambda r: r.protein_id == "protein1")
        ```

        特定の化合物に関連する結果のみを抽出する場合：
        ```python
        compound_results = collection.filter(
            lambda r: r.compound_set_id == "set1" and r.compound_index == 0
        )
        ```

        Args:
            condition: フィルタリング条件

        Returns:
            フィルタリングされた結果のコレクション
        """
        filtered_results = [result for result in self._results if condition(result)]
        filtered_collection = DockingResultCollection(filtered_results)
        return filtered_collection

    def __len__(self) -> int:
        """
        結果の数を取得する。

        Returns:
            結果の数
        """
        return len(self._results)

    def __iter__(self) -> Iterator[DockingResult]:
        """
        結果のイテレータを取得する。

        Returns:
            結果のイテレータ
        """
        return iter(self._results)

    def merge(self, other: "DockingResultCollection") -> "DockingResultCollection":
        """
        別のコレクションとマージして新しいコレクションを作成する。

        Args:
            other: マージする別のコレクション

        Returns:
            マージされた新しいコレクション
        """
        merged = DockingResultCollection()
        for result in itertools.chain(self._results, other._results):
            merged.add(result)
        return merged

    def summarize(self) -> str:
        """
        コレクション内の結果のサマリーを文字列として返す。

        Returns:
            結果のサマリー文字列
        """
        if not self._results:
            return "DockingResultCollection: 結果なし"

        stats = self.get_statistics()
        top_results = self.get_top(5)

        summary = [
            f"DockingResultCollection: {len(self._results)}件の結果",
            f"スコア統計: 最小={stats['min_score']:.3f}, 最大={stats['max_score']:.3f}, 平均={stats['avg_score']:.3f}",
            "上位5件の結果:",
        ]

        for i, result in enumerate(top_results):
            summary.append(
                f"  {i+1}. タンパク質ID={result.protein_id}, 化合物セットID={result.compound_set_id}, "
                f"化合物インデックス={result.compound_index}, スコア={result.docking_score:.3f}"
            )

        return "\n".join(summary)

    def format_top_hits(self, n: int = 5) -> str:
        """
        上位n件の結果をフォーマットした文字列を返す。

        Args:
            n: 表示する上位結果の数

        Returns:
            フォーマットされた上位結果の文字列
        """
        top_hits = self.get_top(n)
        if not top_hits:
            return "上位結果なし"

        lines = [f"上位 {len(top_hits)} 件のドッキング結果:"]
        for i, hit in enumerate(top_hits):
            scores = hit.metadata.get("scores", [])
            score_str = (
                f"全ポーズのスコア: {[float(s) for s in scores[0]]}" if scores else f"スコア: {hit.docking_score:.3f}"
            )
            lines.append(f"  {i+1}. {score_str}")

        return "\n".join(lines)

    def export_to_sdf(self, output_path: Path) -> None:
        """
        結果をSDFファイルにエクスポートする。

        各ドッキング結果の3D構造情報とメタデータ（タンパク質ID、化合物ID、ドッキングスコアなど）を
        含むSDFファイルを生成する。結果はドッキングスコアの昇順でソートされる。

        Args:
            output_path: 出力ファイルのパス
        """
        # 出力ディレクトリが存在しない場合は作成
        output_path.parent.mkdir(parents=True, exist_ok=True)

        # RDKitのSDWriterを使用してSDFファイルを作成
        from rdkit import Chem

        from docking_automation.converters.molecule_converter import MoleculeConverter

        # 空のコレクションの場合でも空のSDFファイルを作成する
        with open(output_path, "w") as f:
            # SDFファイルのヘッダーを書き込む（空のファイルでも有効なSDFファイルとなる）
            f.write("$$$$\n")

        # 結果がある場合のみSDWriterを使用
        if len(self._results) > 0:
            # 処理した結果の数をカウント
            processed_count = 0
            success_count = 0

            with Chem.SDWriter(str(output_path)) as writer:
                # コレクション内の各ドッキング結果を処理
                for result in self._results:
                    processed_count += 1
                    try:
                        # 分子コンバーターを使用して3D構造情報を読み込む
                        converter = MoleculeConverter()
                        mols = converter.pdbqt_to_rdkit(result.result_path)

                        if not mols:
                            continue

                        # 各分子にメタデータを追加して書き込む
                        for mol in mols:
                            # メタデータを追加
                            mol.SetProp("protein_id", result.protein_id)
                            mol.SetProp("compound_set_id", result.compound_set_id)
                            mol.SetProp("compound_index", str(result.compound_index))
                            mol.SetProp("docking_score", str(result.docking_score))

                            # 追加のメタデータを設定
                            for key, value in result.metadata.items():
                                if isinstance(value, (str, int, float, bool)):
                                    mol.SetProp(key, str(value))

                            # SDFファイルに書き込む
                            writer.write(mol)
                            success_count += 1
                    except Exception as e:
                        self._log_warning(f"結果 {result.id} の処理中にエラーが発生しました: {e}")
                        continue

            # 処理結果の概要を表示
            self._log_info(f"SDFエクスポート完了: {success_count}/{processed_count} 件の結果を処理しました。")
        else:
            self._log_warning("エクスポートする結果がありませんでした。")

    def _log_info(self, message: str) -> None:
        """
        情報メッセージをログ出力する（内部使用）。

        Args:
            message: ログメッセージ
        """
        print(f"[DockingResultCollection] {message}")

    def _log_warning(self, message: str) -> None:
        """
        警告メッセージをログ出力する（内部使用）。

        Args:
            message: 警告メッセージ
        """
        print(f"[DockingResultCollection] 警告: {message}")

    @classmethod
    def from_results(cls, results: List[DockingResult]) -> "DockingResultCollection":
        """
        結果のリストからDockingResultCollectionを作成する。

        Args:
            results: DockingResultのリスト

        Returns:
            作成されたDockingResultCollection
        """
        return cls(results)

    @classmethod
    def from_file(cls, file_path: Path) -> "DockingResultCollection":
        """
        ファイルからDockingResultCollectionを作成する。

        Args:
            file_path: 入力ファイルのパス

        Returns:
            作成されたDockingResultCollection
        """
        # 実際のファイル読み込み処理はRDKitなどの外部ライブラリに依存するため、
        # ここではインターフェースのみを定義
        raise NotImplementedError("ファイル読み込み機能を実装する必要があります")

    def get_statistics(self) -> Dict[str, float]:
        """
        コレクション内の結果の統計情報を計算する。

        Returns:
            統計情報を含む辞書（最小値、最大値、平均値、中央値など）
        """
        if not self._results:
            return {"min_score": 0.0, "max_score": 0.0, "avg_score": 0.0, "median_score": 0.0, "count": 0}

        scores = [result.docking_score for result in self._results]
        return {
            "min_score": min(scores),
            "max_score": max(scores),
            "avg_score": sum(scores) / len(scores),
            "median_score": scores[len(scores) // 2],  # 厳密には正確な中央値ではない
            "count": len(scores),
        }
