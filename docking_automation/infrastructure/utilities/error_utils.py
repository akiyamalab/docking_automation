import sys
from contextlib import contextmanager
from io import StringIO
from typing import Generator


@contextmanager
def capture_stderr() -> Generator[StringIO, None, None]:
    """
    標準エラー出力をキャプチャするためのコンテキストマネージャ。

    with句で使用することで、ブロック内での標準エラー出力をキャプチャし、
    ブロック終了時に元の標準エラー出力に戻します。

    Returns:
        キャプチャされた標準エラー出力を含むStringIOオブジェクト

    Examples:
        ```python
        with capture_stderr() as stderr:
            # 何らかの処理（標準エラー出力が発生する可能性がある）
            process_something()

        # キャプチャされた標準エラー出力を取得
        error_output = stderr.getvalue()
        ```
    """
    stderr_buffer = StringIO()
    original_stderr = sys.stderr
    sys.stderr = stderr_buffer
    try:
        yield stderr_buffer
    finally:
        sys.stderr = original_stderr
