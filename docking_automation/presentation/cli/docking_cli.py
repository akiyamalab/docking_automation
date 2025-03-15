#!/usr/bin/env python3
import os
import sys
import argparse
import logging
import json
import configparser
from typing import Dict, Any, Optional

# モックサービスのファクトリをインポート
from docking_automation.presentation.cli.mock_factory import create_mock_services

# 実装が利用可能な場合にロードするためのクラス
# これらはモックモードでは使用されない
# Any型を使用して型チェッカーに対応
from typing import Any, Type, Optional, cast

# サービスクラスの型変数（Any型としてスタブを定義）
OpenBabelFormatConverter: Optional[Any] = None
MGLMoleculePreparationService: Optional[Any] = None
VinaDockingService: Optional[Any] = None

# 実際のサービスが利用可能かチェック
real_services_available = False
try:
    from docking_automation.infrastructure.formats.openbabel_format_converter import OpenBabelFormatConverter
    from docking_automation.infrastructure.tools.mgltools.mgl_molecule_preparation_service import MGLMoleculePreparationService
    from docking_automation.infrastructure.tools.vina.vina_docking_service import VinaDockingService
    from docking_automation.application.workflows.simple_docking_workflow import SimpleDockingWorkflow
    real_services_available = True
except ImportError:
    pass

# SimpleDockingWorkflowは常にインポート
from docking_automation.application.workflows.simple_docking_workflow import SimpleDockingWorkflow

# ロガーの設定
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger("docking_cli")


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """設定ファイルを読み込む
    
    Args:
        config_path: 設定ファイルのパス
        
    Returns:
        設定情報の辞書
    """
    # デフォルト設定
    config = {
        "paths": {
            "obabel": "obabel",
            "prepare_ligand": "prepare_ligand4.py",
            "prepare_receptor": "prepare_receptor4.py",
            "vina": "vina"
        },
        "output_dir": "./output",
        "log_level": "INFO"
    }
    
    # 設定ファイルがあれば読み込む
    if config_path and os.path.exists(config_path):
        try:
            parser = configparser.ConfigParser()
            parser.read(config_path)
            
            # パスセクション
            if "PATHS" in parser:
                # 辞書に明示的なキャストを使用して型を指定
                paths_dict = cast(Dict[str, str], config["paths"])
                paths_section = cast(Dict[str, str], parser["PATHS"])
                for key in paths_section:
                    paths_dict[key] = paths_section[key]
            
            # その他の設定
            if "GENERAL" in parser:
                if "output_dir" in parser["GENERAL"]:
                    config["output_dir"] = parser["GENERAL"]["output_dir"]
                if "log_level" in parser["GENERAL"]:
                    config["log_level"] = parser["GENERAL"]["log_level"]
        
        except Exception as e:
            logger.error(f"Error loading config file: {e}")
    
    return config

def setup_services(config: Dict[str, Any], mock_mode: bool = False) -> SimpleDockingWorkflow:
    """サービスをセットアップする
    
    Args:
        config: 設定情報
        mock_mode: モックモードで実行するかどうか
        
    Returns:
        セットアップされたワークフロー
    """
    # モードの強制設定をチェック
    use_mock = os.environ.get("USE_MOCK_SERVICES") == "1" or mock_mode
    
    # 実際のサービスが利用可能かつモックモードが指定されていない場合
    if not use_mock and real_services_available and OpenBabelFormatConverter and MGLMoleculePreparationService and VinaDockingService:
        logger.info("Using real services")
        try:
            # パスの取得
            obabel_path = config["paths"]["obabel"]
            prepare_ligand_path = config["paths"]["prepare_ligand"]
            prepare_receptor_path = config["paths"]["prepare_receptor"]
            vina_path = config["paths"]["vina"]
            
            # サービスの作成（実クラスを使用）
            fc = OpenBabelFormatConverter  # 型チェックのための一時変数
            mps = MGLMoleculePreparationService  # 型チェックのための一時変数
            ds = VinaDockingService  # 型チェックのための一時変数
            
            format_converter = fc(obabel_path=obabel_path)
            molecule_preparation_service = mps(
                prepare_ligand_path=prepare_ligand_path,
                prepare_receptor_path=prepare_receptor_path,
                obabel_path=obabel_path
            )
            docking_service = ds(vina_path=vina_path)
            
            # ワークフローの作成
            workflow = SimpleDockingWorkflow(
                molecule_preparation_service=molecule_preparation_service,
                docking_service=docking_service
            )
            
            return workflow
        except Exception as e:
            logger.warning(f"Failed to initialize real services: {e}")
            logger.warning("Falling back to mock services")
    
    # モックサービスを使用
    logger.info("Using mock services")
    services = create_mock_services(config)
    workflow = services["workflow"]
    assert isinstance(workflow, SimpleDockingWorkflow)
    return workflow


def parse_args() -> argparse.Namespace:
    """コマンドライン引数を解析する"""
    parser = argparse.ArgumentParser(description="Docking Automation Tool")
    
    # サブコマンドの追加
    subparsers = parser.add_subparsers(dest="command", help="Command to execute")
    
    # dockコマンド
    dock_parser = subparsers.add_parser("dock", help="Perform docking calculation")
    dock_parser.add_argument("-l", "--ligand", required=True, help="Ligand file path")
    dock_parser.add_argument("-r", "--receptor", required=True, help="Receptor file path")
    dock_parser.add_argument("-c", "--config", help="Configuration file path")
    dock_parser.add_argument("-o", "--output-dir", help="Output directory path")
    dock_parser.add_argument("-p", "--parameters", help="Docking parameters in JSON format")
    dock_parser.add_argument("-m", "--mock", action="store_true", help="Use mock implementation (no external tools required)")
    
    # バージョン表示
    parser.add_argument("-v", "--version", action="store_true", help="Show version")
    
    return parser.parse_args()


def main() -> int:
    """メイン関数"""
    try:
        # コマンドライン引数の解析
        args = parse_args()
        
        # バージョン表示
        if args.version:
            print("Docking Automation Tool v0.1.0")
            return 0
        
        # コマンドがない場合はヘルプを表示
        if not args.command:
            print("Please specify a command. Use -h for help.")
            return 1
        
        # 設定の読み込み
        config = load_config(args.config if hasattr(args, "config") else None)
        
        # ログレベルの設定
        log_level = getattr(logging, config["log_level"], logging.INFO)
        logging.getLogger().setLevel(log_level)
        
        # 出力ディレクトリの設定
        output_dir = args.output_dir if hasattr(args, "output_dir") and args.output_dir else config["output_dir"]
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # サービスのセットアップ（mock引数があれば渡す）
        mock_mode = False
        if hasattr(args, "mock") and args.mock:
            mock_mode = True
            
        workflow = setup_services(config, mock_mode)
        
        # コマンドの実行
        if args.command == "dock":
            # ドッキングパラメータの解析
            parameters = None
            if hasattr(args, "parameters") and args.parameters:
                try:
                    parameters = json.loads(args.parameters)
                except json.JSONDecodeError as e:
                    logger.error(f"Error parsing parameters: {e}")
                    return 1
            
            # ドッキング計算の実行
            logger.info(f"Starting docking calculation")
            logger.info(f"Ligand: {args.ligand}")
            logger.info(f"Receptor: {args.receptor}")
            
            result = workflow.execute(
                ligand_path=args.ligand,
                receptor_path=args.receptor,
                config=parameters,
                output_dir=output_dir
            )
            
            logger.info(f"Docking calculation completed")
            best_score = result.get_best_score()
            if best_score:
                logger.info(f"Best score: {best_score.value} {best_score.unit if best_score.unit else ''}")
            logger.info(f"Number of poses: {len(result.poses)}")
            logger.info(f"Results saved to: {output_dir}")
            
        return 0
        
    except Exception as e:
        logger.error(f"Error: {e}")
        return 1


if __name__ == "__main__":
    sys.exit(main())