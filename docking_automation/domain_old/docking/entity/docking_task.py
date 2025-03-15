import uuid
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Union
from enum import Enum, auto

from docking_automation.domain.docking.entity.ligand import Ligand
from docking_automation.domain.docking.entity.receptor import Receptor
from docking_automation.domain.docking.value_object.docking_configuration import DockingConfiguration


class TaskStatus(Enum):
    """ドッキングタスクのステータスを表す列挙型"""
    PENDING = auto()      # 保留中
    PREPARING = auto()    # 準備中
    READY = auto()        # 準備完了、実行待ち
    RUNNING = auto()      # 実行中
    COMPLETED = auto()    # 完了
    FAILED = auto()       # 失敗
    CANCELLED = auto()    # キャンセル


@dataclass
class DockingTask:
    """ドッキング計算タスクを表すエンティティ
    
    リガンド、レセプター、ドッキング設定を含み、ドッキング計算のコンテキストを表します。
    """
    
    id: str
    ligand: Ligand
    receptor: Union[Receptor, List[Receptor]]  # 単一または複数のレセプターをサポート
    configuration: DockingConfiguration
    status: TaskStatus = TaskStatus.PENDING
    error_message: Optional[str] = None
    metadata: Dict[str, Any] = field(default_factory=dict)
    created_at: Optional[float] = None  # Unix timestamp
    updated_at: Optional[float] = None  # Unix timestamp
    timeout_seconds: int = 3600  # デフォルトタイムアウト：1時間
    molecular_weight: Optional[float] = None  # 分子量（フィルタリング用）
    
    def __post_init__(self) -> None:
        """初期化後の処理"""
        import time
        
        # IDが指定されていない場合はUUIDを生成
        if self.id == "":
            self.id = uuid.uuid4().hex
        
        # 作成日時が指定されていない場合は現在時刻を設定
        if self.created_at is None:
            self.created_at = time.time()
            self.updated_at = self.created_at
    
    def update_status(self, status: TaskStatus, error_message: Optional[str] = None) -> None:
        """タスクのステータスを更新
        
        Args:
            status: 新しいステータス
            error_message: エラーメッセージ（失敗時）
        """
        import time
        
        self.status = status
        self.updated_at = time.time()
        
        if error_message:
            self.error_message = error_message
    def get_receptors(self) -> List[Receptor]:
        """単一または複数のレセプターをリストとして取得"""
        if isinstance(self.receptor, list):
            return self.receptor
        return [self.receptor]
    
    def is_multi_receptor(self) -> bool:
        """複数のレセプターを持つかどうか"""
        return isinstance(self.receptor, list) and len(self.receptor) > 1
    
    def get_primary_receptor(self) -> Receptor:
        """主要なレセプターを取得（複数ある場合は最初のもの）"""
        if isinstance(self.receptor, list):
            if self.receptor:
                return self.receptor[0]
            raise ValueError("No receptors available in this task")
        return self.receptor
        
    def is_ready(self) -> bool:
        """タスクが実行可能かどうかを確認
        
        リガンドとレセプターが適切に準備されているか確認します。
        """
        # リガンドが準備されているか確認
        if not self.ligand.is_prepared():
            return False
            
        # すべてのレセプターが準備されているか確認
        receptors = self.get_receptors()
        for receptor in receptors:
            if not receptor.is_prepared():
                return False
                
        # 設定が有効か確認
        if not self.configuration.validate():
            return False
            
        # 分子量フィルターがある場合は確認
        if self.molecular_weight is not None and self.molecular_weight >= 700:
            return False
            
        return True
    
    def prepare_for_execution(self) -> bool:
        """タスクを実行するための準備を行う
        
        ステータスをPREPARINGに更新し、リガンドとレセプターが準備されていることを確認します。
        準備が完了したらステータスをREADYに更新します。
        
        Returns:
            準備が成功したかどうか
        """
        # 既にREADYまたはそれ以降のステータスであれば何もしない
        if self.status in [TaskStatus.READY, TaskStatus.RUNNING, TaskStatus.COMPLETED]:
            return True
        
        # 準備中にする
        self.update_status(TaskStatus.PREPARING)
        
        # リガンドとレセプターが準備されているかチェック
        if not self.is_ready():
            self.update_status(TaskStatus.FAILED, "Ligand or receptor is not properly prepared")
            return False
        
        # 準備完了
        self.update_status(TaskStatus.READY)
        return True
    
    def mark_as_running(self) -> None:
        """タスクを実行中にマーク"""
        self.update_status(TaskStatus.RUNNING)
    
    def mark_as_completed(self) -> None:
        """タスクを完了にマーク"""
        self.update_status(TaskStatus.COMPLETED)
    
    def mark_as_failed(self, error_message: str) -> None:
        """タスクを失敗にマーク"""
        self.update_status(TaskStatus.FAILED, error_message)
    
    def mark_as_cancelled(self) -> None:
        """タスクをキャンセルにマーク"""
        self.update_status(TaskStatus.CANCELLED)
    
    def is_active(self) -> bool:
        """タスクがアクティブかどうかを確認"""
        return self.status in [TaskStatus.PREPARING, TaskStatus.READY, TaskStatus.RUNNING]
    
    def is_finished(self) -> bool:
        """タスクが完了したかどうかを確認（成功、失敗、キャンセルのいずれか）"""
        return self.status in [TaskStatus.COMPLETED, TaskStatus.FAILED, TaskStatus.CANCELLED]
    
    def get_time_elapsed(self) -> Optional[float]:
        """タスクの経過時間を取得（秒）"""
        import time
        
        if not self.created_at:
            return None
        
        if self.is_finished() and self.updated_at:
            return self.updated_at - self.created_at
        
        return time.time() - self.created_at
    
    def is_timed_out(self) -> bool:
        """タスクがタイムアウトしたかどうかを確認"""
        elapsed = self.get_time_elapsed()
        return elapsed is not None and elapsed > self.timeout_seconds
    
    def __str__(self) -> str:
        """文字列表現"""
        receptor_info = ""
        if isinstance(self.receptor, list):
            receptor_count = len(self.receptor)
            if receptor_count > 0:
                primary_receptor = self.receptor[0]
                primary_name = primary_receptor.name if primary_receptor.name is not None else "Unnamed"
                receptor_info = f"{primary_name} (+{receptor_count-1} more)" if receptor_count > 1 else primary_name
            else:
                receptor_info = "None"
        else:
            receptor_info = self.receptor.name if self.receptor.name is not None else "Unnamed"
            
        return f"DockingTask(id={self.id}, ligand={self.ligand.name}, receptor={receptor_info}, status={self.status.name})"