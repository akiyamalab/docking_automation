from abc import ABC, abstractmethod
from typing import Optional, List, Dict, Any

from docking_automation.docking.entity.docking_task import DockingTask
from docking_automation.docking.entity.docking_result import DockingResult
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration


class DockingService(ABC):
    """ドッキング計算を実行するサービスのインターフェース
    
    このサービスは、リガンドとレセプターを受け取り、指定された設定でドッキング計算を実行し、
    結果を返す責務を持ちます。具体的な実装はツールごとに異なります。
    """
    
    @abstractmethod
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
        pass
    
    @abstractmethod
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
        pass
    
    @abstractmethod
    def validate_task(self, task: DockingTask) -> bool:
        """タスクの妥当性を検証
        
        Args:
            task: 検証するドッキングタスク
            
        Returns:
            タスクが有効かどうか
        """
        pass
    
    @abstractmethod
    def get_supported_parameters(self) -> List[str]:
        """このサービスがサポートするパラメータの一覧を取得
        
        Returns:
            サポートされているパラメータ名のリスト
        """
        pass
    
    @abstractmethod
    def get_default_configuration(self) -> DockingConfiguration:
        """このサービスのデフォルト設定を取得
        
        Returns:
            デフォルトのドッキング設定
        """
        pass
    
    @abstractmethod
    def cancel_task(self, task: DockingTask) -> bool:
        """実行中のタスクをキャンセル
        
        Args:
            task: キャンセルするタスク
            
        Returns:
            キャンセルが成功したかどうか
        """
        pass
    
    @abstractmethod
    def get_result_for_task(self, task_id: str) -> Optional[DockingResult]:
        """タスクIDに対応する結果を取得（キャッシュから）
        
        Args:
            task_id: 取得するタスクのID
            
        Returns:
            タスクの結果（見つからない場合はNone）
        """
        pass