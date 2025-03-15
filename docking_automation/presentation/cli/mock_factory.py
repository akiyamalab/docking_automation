"""モックサービスのファクトリモジュール

このモジュールは、外部ツールがインストールされていない環境でも動作する
モックサービスを提供するファクトリを定義します。
"""

import logging
from typing import Dict, Any, Optional

from docking_automation.molecule.service.molecule_preparation_service import MoleculePreparationService
from docking_automation.docking.service.docking_service import DockingService
from docking_automation.infrastructure.tools.mock.mock_molecule_preparation_service import MockMoleculePreparationService
from docking_automation.infrastructure.tools.mock.mock_docking_service import MockDockingService
from docking_automation.application.workflows.simple_docking_workflow import SimpleDockingWorkflow

# ロガーの設定
logger = logging.getLogger(__name__)


def create_mock_services(config: Optional[Dict[str, Any]] = None) -> Dict[str, Any]:
    """モックサービスを作成する
    
    Args:
        config: 設定情報（使用されません）
        
    Returns:
        サービスの辞書
    """
    logger.info("Creating mock services")
    
    # モックサービスのインスタンス化
    molecule_preparation_service = MockMoleculePreparationService()
    docking_service = MockDockingService()
    
    # ワークフローの作成
    workflow = SimpleDockingWorkflow(
        molecule_preparation_service=molecule_preparation_service,
        docking_service=docking_service
    )
    
    return {
        "molecule_preparation_service": molecule_preparation_service,
        "docking_service": docking_service,
        "workflow": workflow
    }