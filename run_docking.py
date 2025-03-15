#!/usr/bin/env python3
import os
import sys
import logging

# 現在のファイルの絶対パスを取得
current_dir = os.path.dirname(os.path.abspath(__file__))

# パッケージのルートをPYTHONPATHに追加（明示的にパスを指定）
sys.path.insert(0, current_dir)

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("docking_script")
def main():
    try:
        # パス情報を出力（デバッグ用）
        logger.info(f"Current directory: {current_dir}")
        logger.info(f"Python path: {sys.path}")
        
        # モック実装を強制的に使用するために環境変数を設定
        os.environ["USE_MOCK_SERVICES"] = "1"
        
        # CLIモジュールからメイン関数をインポート
        from docking_automation.presentation.cli.docking_cli import main as cli_main
        
        # コマンドライン引数を修正して--mockオプションを追加
        if len(sys.argv) > 1 and sys.argv[1] == "dock":
            if "--mock" not in sys.argv and "-m" not in sys.argv:
                sys.argv.append("--mock")
        
        # メイン関数を実行
        return cli_main()
        return cli_main()
    
    except ImportError as e:
        logger.error(f"Import error: {e}")
        logger.error("Failed to import docking modules. Make sure the package is installed correctly.")
        return 1
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())