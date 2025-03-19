from typing import Dict, List, Type, Callable, Set
from docking_automation.domain.domain_event import DomainEvent

class DomainEventPublisher:
    """
    ドメインイベントの発行と購読を管理するクラス。
    
    シングルトンパターンを使用して、アプリケーション全体で一貫したイベント管理を提供する。
    """
    _subscribers: Dict[Type[DomainEvent], List[Callable[[DomainEvent], None]]] = {}
    
    @classmethod
    def subscribe(cls, event_type: Type[DomainEvent], subscriber: Callable[[DomainEvent], None]) -> None:
        """
        イベントタイプに対する購読者を登録する。
        
        Args:
            event_type: 購読するイベントの型
            subscriber: イベント発行時に呼び出されるコールバック関数
        """
        if event_type not in cls._subscribers:
            cls._subscribers[event_type] = []
        cls._subscribers[event_type].append(subscriber)
    
    @classmethod
    def unsubscribe(cls, event_type: Type[DomainEvent], subscriber: Callable[[DomainEvent], None]) -> None:
        """
        イベントタイプに対する購読者の登録を解除する。
        
        Args:
            event_type: 購読解除するイベントの型
            subscriber: 登録解除するコールバック関数
        """
        if event_type in cls._subscribers and subscriber in cls._subscribers[event_type]:
            cls._subscribers[event_type].remove(subscriber)
    
    @classmethod
    def publish(cls, event: DomainEvent) -> None:
        """
        イベントを発行し、購読者に通知する。
        
        Args:
            event: 発行するイベント
        """
        event_type = type(event)
        
        # 特定のイベント型の購読者に通知
        if event_type in cls._subscribers:
            for subscriber in cls._subscribers[event_type]:
                subscriber(event)
        
        # DomainEventの購読者にも通知（全てのイベントを購読）
        if DomainEvent in cls._subscribers and event_type != DomainEvent:
            for subscriber in cls._subscribers[DomainEvent]:
                subscriber(event)
    
    @classmethod
    def reset(cls) -> None:
        """
        全ての購読者登録をリセットする。
        主にテスト用途で使用する。
        """
        cls._subscribers.clear()