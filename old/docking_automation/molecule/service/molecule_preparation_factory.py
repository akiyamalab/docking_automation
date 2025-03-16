"""分子準備サービスのファクトリークラス

このモジュールは、各ドッキングツールに適した分子準備サービスを生成するファクトリークラスを提供します。
"""

from typing import Dict, Type

from docking_automation.molecule.service.molecule_preparation_service import MoleculePreparationService
from docking_automation.molecule.service.vina_preparation_service import VinaPreparationService


class MoleculePreparationFactory:
    """分子準備サービスのファクトリークラス
    
    各ドッキングツールに適した分子準備サービスのインスタンスを生成します。
    """
    
    # ドッキングツール名と対応する準備サービスクラスのマッピング
    _service_classes: Dict[str, Type[MoleculePreparationService]] = {
        "vina": VinaPreparationService,
        # 他のサービスは必要に応じて追加
    }
    
    @classmethod
    def create_for_tool(cls, tool_name: str) -> MoleculePreparationService:
        """指定されたドッキングツール用の分子準備サービスを生成する
        
        Args:
            tool_name: ドッキングツール名 ("vina", "gnina", "restretto")
            
        Returns:
            適切な分子準備サービスのインスタンス
            
        Raises:
            ValueError: サポートされていないドッキングツール名が指定された場合
        """
        tool_name = tool_name.lower()
        
        if tool_name not in cls._service_classes:
            supported_tools = ", ".join(cls._service_classes.keys())
            raise ValueError(
                f"未サポートのドッキングツール: {tool_name}. "
                f"サポートされているツール: {supported_tools}"
            )
        
        service_class = cls._service_classes[tool_name]
        return service_class()
    
    @classmethod
    def get_supported_tools(cls) -> Dict[str, str]:
        """サポートされているドッキングツール名の一覧を取得する
        
        Returns:
            ドッキングツール名とその説明のマッピング
        """
        return {
            "vina": "AutoDock Vina: 高速で精度の高いオープンソースのドッキングツール",
            # 他のツールは必要に応じて追加
        }