#!/usr/bin/env python3
import os
import sys
import logging

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("docking_script")

def main() -> int:
    try:
        # 直接CLIモジュールからmain関数をインポート
        from docking_automation.presentation.cli.docking_cli import main as cli_main
        
        # メイン関数を実行
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