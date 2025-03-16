#!/usr/bin/env python3
"""
DockingConfigurationクラスのテスト

このモジュールはドッキング計算の設定を表すDockingConfigurationクラスをテストします。
DockingConfigurationはグリッドボックスやパラメータなど、ドッキング計算に必要な設定情報を保持します。
"""

import pytest
from dataclasses import FrozenInstanceError
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.docking.value_object.docking_parameter import DockingParameters, DockingParameter, ParameterType


class TestDockingConfiguration:
    """DockingConfigurationクラスのテストケース"""
    
    @pytest.fixture
    def sample_grid_box(self):
        """テスト用のグリッドボックス"""
        return GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
    
    @pytest.fixture
    def sample_parameters(self):
        """テスト用のドッキングパラメータ"""
        params = DockingParameters()
        params.add(DockingParameter(
            name="exhaustiveness",
            value=8,
            parameter_type=ParameterType.INTEGER
        ))
        params.add(DockingParameter(
            name="num_modes",
            value=9,
            parameter_type=ParameterType.INTEGER
        ))
        params.add(DockingParameter(
            name="energy_range",
            value=3.0,
            parameter_type=ParameterType.FLOAT
        ))
        return params
    
    def test_initialization(self, sample_grid_box, sample_parameters):
        """初期化と基本プロパティのテスト"""
        # 最小限の引数での初期化
        config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters
        )
        
        assert config.grid_box is sample_grid_box
        assert config.parameters is sample_parameters
        assert config.name is None
        assert config.description is None
        assert config.metadata == {}
        
        # 全オプション付きでの初期化
        metadata = {"tool": "vina", "version": "1.2.0"}
        config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters,
            name="Test Configuration",
            description="A test configuration for unit testing",
            metadata=metadata
        )
        
        assert config.grid_box is sample_grid_box
        assert config.parameters is sample_parameters
        assert config.name == "Test Configuration"
        assert config.description == "A test configuration for unit testing"
        assert config.metadata == metadata
    
    def test_immutability(self, sample_grid_box, sample_parameters):
        """イミュータブル性のテスト"""
        config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters,
            name="Test Configuration"
        )
        
        # 直接プロパティを変更するとエラーになる
        with pytest.raises(FrozenInstanceError):
            config.name = "Modified Configuration"
        
        with pytest.raises(FrozenInstanceError):
            config.grid_box = GridBox(
                center_x=0.0,
                center_y=0.0,
                center_z=0.0,
                size_x=10.0,
                size_y=10.0,
                size_z=10.0
            )
        
        with pytest.raises(FrozenInstanceError):
            config.parameters = DockingParameters()
        
        with pytest.raises(FrozenInstanceError):
            config.metadata = {"new": "value"}
    
    def test_validate(self, sample_grid_box, sample_parameters):
        """validate()メソッドのテスト"""
        # 有効な設定
        valid_config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters
        )
        assert valid_config.validate() is True
        
        # 無効なグリッドボックスを持つ設定をテスト
        # グリッドボックスvalidationが失敗するケース（サイズに負の値）
        with pytest.raises(ValueError):
            # このグリッドボックス生成自体がエラーになるはず
            invalid_grid = GridBox(
                center_x=10.0,
                center_y=20.0,
                center_z=30.0,
                size_x=-15.0,  # 負の値
                size_y=25.0,
                size_z=35.0
            )
        
        # validateメソッドのパラメータ検証部分のテスト
        # monkeypatchを使ってDockingParametersのvalidateメソッドを一時的に差し替え
        invalid_params = sample_parameters
        original_validate = invalid_params.validate
        try:
            # validateメソッドをオーバーライド
            invalid_params.validate = lambda: False
            
            # 無効なパラメータを持つ設定
            invalid_config = DockingConfiguration(
                grid_box=sample_grid_box,
                parameters=invalid_params
            )
            assert invalid_config.validate() is False
        finally:
            # 元に戻す
            invalid_params.validate = original_validate
    
    def test_with_parameters(self, sample_grid_box, sample_parameters):
        """with_parameters()メソッドのテスト"""
        original_config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters,
            name="Original Configuration"
        )
        
        # 新しいパラメータ
        new_params = DockingParameters()
        new_params.add(DockingParameter(
            name="new_param",
            value="new_value",
            parameter_type=ParameterType.STRING
        ))
        
        # 新しいパラメータで設定を更新
        updated_config = original_config.with_parameters(new_params)
        
        # 元の設定は変更されていないことを確認
        assert original_config.parameters is sample_parameters
        assert "new_param" not in [param.name for param in original_config.parameters.parameters.values()]
        
        # 更新された設定は新しいインスタンスであることを確認
        assert updated_config is not original_config
        
        # グリッドボックスは同じであることを確認
        assert updated_config.grid_box is original_config.grid_box
        
        # その他のプロパティも同じであることを確認
        assert updated_config.name == original_config.name
        assert updated_config.description == original_config.description
        assert updated_config.metadata == original_config.metadata
        
        # パラメータがマージされていることを確認
        assert updated_config.parameters is not original_config.parameters
        assert updated_config.parameters is not new_params
        assert updated_config.get_parameter("new_param") == "new_value"
        
        # 元のパラメータも保持されていることを確認
        assert updated_config.get_parameter("exhaustiveness") == 8
        assert updated_config.get_parameter("num_modes") == 9
        assert updated_config.get_parameter("energy_range") == 3.0
    
    def test_with_grid_box(self, sample_grid_box, sample_parameters):
        """with_grid_box()メソッドのテスト"""
        original_config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters,
            name="Original Configuration"
        )
        
        # 新しいグリッドボックス
        new_grid_box = GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=10.0,
            size_y=10.0,
            size_z=10.0
        )
        
        # 新しいグリッドボックスで設定を更新
        updated_config = original_config.with_grid_box(new_grid_box)
        
        # 元の設定は変更されていないことを確認
        assert original_config.grid_box is sample_grid_box
        
        # 更新された設定は新しいインスタンスであることを確認
        assert updated_config is not original_config
        
        # グリッドボックスが更新されていることを確認
        assert updated_config.grid_box is new_grid_box
        
        # パラメータは同じであることを確認
        assert updated_config.parameters is original_config.parameters
        
        # その他のプロパティも同じであることを確認
        assert updated_config.name == original_config.name
        assert updated_config.description == original_config.description
        assert updated_config.metadata == original_config.metadata
    
    def test_get_parameter(self, sample_grid_box, sample_parameters):
        """get_parameter()メソッドのテスト"""
        config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters
        )
        
        # 存在するパラメータの取得
        assert config.get_parameter("exhaustiveness") == 8
        assert config.get_parameter("num_modes") == 9
        assert config.get_parameter("energy_range") == 3.0
        
        # 存在しないパラメータの取得（デフォルトなし）
        assert config.get_parameter("non_existent") is None
        
        # 存在しないパラメータの取得（デフォルトあり）
        assert config.get_parameter("non_existent", default=42) == 42
    
    def test_to_dict(self, sample_grid_box, sample_parameters):
        """to_dict()メソッドのテスト"""
        # 最小限のプロパティでの変換
        min_config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters
        )
        
        min_dict = min_config.to_dict()
        assert "grid_box" in min_dict
        assert "parameters" in min_dict
        assert "name" not in min_dict
        assert "description" not in min_dict
        assert "metadata" not in min_dict
        
        # グリッドボックスとパラメータが正しく変換されていることを確認
        assert min_dict["grid_box"]["center_x"] == 10.0
        assert min_dict["grid_box"]["center_y"] == 20.0
        assert min_dict["grid_box"]["center_z"] == 30.0
        assert min_dict["parameters"]["exhaustiveness"] == 8
        assert min_dict["parameters"]["num_modes"] == 9
        assert min_dict["parameters"]["energy_range"] == 3.0
        
        # 全プロパティでの変換
        full_config = DockingConfiguration(
            grid_box=sample_grid_box,
            parameters=sample_parameters,
            name="Full Configuration",
            description="A configuration with all properties",
            metadata={"key": "value"}
        )
        
        full_dict = full_config.to_dict()
        assert "grid_box" in full_dict
        assert "parameters" in full_dict
        assert "name" in full_dict
        assert "description" in full_dict
        assert "metadata" in full_dict
        assert full_dict["name"] == "Full Configuration"
        assert full_dict["description"] == "A configuration with all properties"
        assert full_dict["metadata"] == {"key": "value"}
    
    def test_from_dict(self):
        """from_dict()メソッドのテスト"""
        # 入力辞書
        input_dict = {
            "grid_box": {
                "center_x": 10.0,
                "center_y": 20.0,
                "center_z": 30.0,
                "size_x": 15.0,
                "size_y": 25.0,
                "size_z": 35.0
            },
            "parameters": {
                "exhaustiveness": 8,
                "num_modes": 9,
                "energy_range": 3.0
            },
            "name": "Test Configuration",
            "description": "A test configuration",
            "metadata": {"key": "value"}
        }
        
        # 辞書から設定の作成
        config = DockingConfiguration.from_dict(input_dict)
        
        # グリッドボックスの検証
        assert config.grid_box.center_x == 10.0
        assert config.grid_box.center_y == 20.0
        assert config.grid_box.center_z == 30.0
        assert config.grid_box.size_x == 15.0
        assert config.grid_box.size_y == 25.0
        assert config.grid_box.size_z == 35.0
        
        # パラメータの検証
        assert config.get_parameter("exhaustiveness") == 8
        assert config.get_parameter("num_modes") == 9
        assert config.get_parameter("energy_range") == 3.0
        
        # その他のプロパティの検証
        assert config.name == "Test Configuration"
        assert config.description == "A test configuration"
        assert config.metadata == {"key": "value"}
        
        # 最小限の辞書からの作成
        min_dict = {
            "grid_box": {
                "center_x": 0.0,
                "center_y": 0.0,
                "center_z": 0.0,
                "size_x": 10.0,
                "size_y": 10.0,
                "size_z": 10.0
            },
            "parameters": {}
        }
        
        min_config = DockingConfiguration.from_dict(min_dict)
        assert min_config.name is None
        assert min_config.description is None
        assert min_config.metadata == {}
    
    def test_create_default(self):
        """create_default()メソッドのテスト"""
        # デフォルト設定の作成
        config = DockingConfiguration.create_default(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0
        )
        
        # グリッドボックスの検証
        assert config.grid_box.center_x == 10.0
        assert config.grid_box.center_y == 20.0
        assert config.grid_box.center_z == 30.0
        assert config.grid_box.size_x == 20.0  # デフォルト値
        assert config.grid_box.size_y == 20.0  # デフォルト値
        assert config.grid_box.size_z == 20.0  # デフォルト値
        
        # パラメータの検証（AutoDock Vinaのデフォルトパラメータ）
        assert config.get_parameter("exhaustiveness") == 8
        assert config.get_parameter("num_modes") == 9
        assert config.get_parameter("energy_range") == 3
        assert config.get_parameter("seed") == 0
        assert config.get_parameter("cpu") == 1
        
        # 名前の検証
        assert config.name == "Default Configuration"
        
        # 他のプロパティはデフォルト値
        assert config.description is None
        assert config.metadata == {}