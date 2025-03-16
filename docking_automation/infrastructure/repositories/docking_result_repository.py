from docking_automation.docking.docking_result import DockingResult


class DockingResultRepository:
    """ドッキング計算結果を保持するリポジトリクラス。"""
    def save(self, result: DockingResult):
        """ドッキング計算結果を保存する。"""
        ...
    def load(self, result_id: str) -> DockingResult:
        """ドッキング計算結果を読み込む。"""
        ...
