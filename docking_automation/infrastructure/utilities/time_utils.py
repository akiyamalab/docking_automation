"""
時間計測に関するユーティリティモジュール。

このモジュールは、関数やメソッドの実行時間を計測するためのデコレータなどを提供します。
"""

import functools
import time

# 循環参照を避けるため、型チェックのみに使用する型インポート
from typing import TYPE_CHECKING, Any, Callable, Optional, TypeVar, cast

if TYPE_CHECKING:
    from docking_automation.docking.docking_result_collection import (
        DockingResultCollection,
    )

T = TypeVar("T")


def measure_execution_time(func: Callable[..., T]) -> Callable[..., T]:
    """
    関数の実行時間を計測し、戻り値のDockingResultCollectionオブジェクトの
    execution_timeフィールドに設定するデコレータ。

    Args:
        func: 計測対象の関数（DockingResultCollectionを返すことを想定）

    Returns:
        実行時間が設定されたDockingResultCollectionを返す新しい関数

    Examples:
        ```python
        @measure_execution_time
        def run_docking(...) -> DockingResultCollection:
            # 処理
            return results
        ```
    """

    @functools.wraps(func)
    def wrapper(*args: Any, **kwargs: Any) -> T:
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        execution_time = end_time - start_time

        # 戻り値がDockingResultCollectionインスタンスであることを確認
        # 実行時にインポートして循環参照を回避
        from docking_automation.docking.docking_result_collection import (
            DockingResultCollection,
        )

        if isinstance(result, DockingResultCollection):
            result.execution_time = execution_time
            print(f"[時間計測] {func.__name__} の実行時間: {execution_time:.4f} 秒")
        else:
            # 期待しない型が返された場合の警告
            print(f"警告: 関数 {func.__name__} はDockingResultCollectionを返しませんでした。実行時間は設定されません。")

        return result

    return wrapper
