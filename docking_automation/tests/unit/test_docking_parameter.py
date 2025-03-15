#!/usr/bin/env python3
"""
DockingParameterクラスのテスト

このモジュールはドッキング計算のパラメータを表すDockingParameterクラスと
パラメータコレクションを表すDockingParametersクラスをテストします。
"""

import pytest
from docking_automation.docking.value_object.docking_parameter import (
    DockingParameter, DockingParameters, ParameterType
)


class TestParameterType:
    """ParameterType列挙型のテストケース"""
    
    def test_from_value(self):
        """値からパラメータ型を推測するテスト"""
        # 整数
        assert ParameterType.from_value(10) == ParameterType.INTEGER
        assert ParameterType.from_value(-5) == ParameterType.INTEGER
        
        # 浮動小数点数
        assert ParameterType.from_value(3.14) == ParameterType.FLOAT
        assert ParameterType.from_value(-0.5) == ParameterType.FLOAT
        
        # 論理値
        assert ParameterType.from_value(True) == ParameterType.BOOLEAN
        assert ParameterType.from_value(False) == ParameterType.BOOLEAN
        
        # 文字列
        assert ParameterType.from_value("text") == ParameterType.STRING
        assert ParameterType.from_value("") == ParameterType.STRING


class TestDockingParameter:
    """DockingParameterクラスのテストケース"""
    
    def test_initialization(self):
        """初期化と基本プロパティのテスト"""
        # 整数パラメータ
        int_param = DockingParameter(
            name="exhaustiveness",
            value=8,
            parameter_type=ParameterType.INTEGER,
            description="Search thoroughness"
        )
        
        assert int_param.name == "exhaustiveness"
        assert int_param.value == 8
        assert int_param.parameter_type == ParameterType.INTEGER
        assert int_param.description == "Search thoroughness"
        assert int_param.min_value is None
        assert int_param.max_value is None
        
        # 浮動小数点パラメータ（範囲制限付き）
        float_param = DockingParameter(
            name="energy_range",
            value=4.0,
            parameter_type=ParameterType.FLOAT,
            min_value=0.0,
            max_value=10.0
        )
        
        assert float_param.name == "energy_range"
        assert float_param.value == 4.0
        assert float_param.parameter_type == ParameterType.FLOAT
        assert float_param.min_value == 0.0
        assert float_param.max_value == 10.0
    
    def test_type_validation(self):
        """パラメータ型のバリデーション"""
        # 有効な型
        DockingParameter(name="int_param", value=5, parameter_type=ParameterType.INTEGER)
        DockingParameter(name="float_param", value=3.14, parameter_type=ParameterType.FLOAT)
        DockingParameter(name="bool_param", value=True, parameter_type=ParameterType.BOOLEAN)
        DockingParameter(name="str_param", value="text", parameter_type=ParameterType.STRING)
        
        # 無効な型（整数に文字列）
        with pytest.raises(ValueError) as excinfo:
            DockingParameter(name="invalid_int", value="text", parameter_type=ParameterType.INTEGER)
        assert "must be an integer" in str(excinfo.value)
        
        # 無効な型（浮動小数点に文字列）
        with pytest.raises(ValueError) as excinfo:
            DockingParameter(name="invalid_float", value="text", parameter_type=ParameterType.FLOAT)
        assert "must be a float" in str(excinfo.value)
        
        # 無効な型（論理値に整数）
        with pytest.raises(ValueError) as excinfo:
            DockingParameter(name="invalid_bool", value=1, parameter_type=ParameterType.BOOLEAN)
        assert "must be a boolean" in str(excinfo.value)
    
    def test_range_validation(self):
        """値の範囲検証"""
        # 範囲内の値
        DockingParameter(
            name="valid_range",
            value=5,
            parameter_type=ParameterType.INTEGER,
            min_value=0,
            max_value=10
        )
        
        # 最小値未満
        with pytest.raises(ValueError) as excinfo:
            DockingParameter(
                name="below_min",
                value=-1,
                parameter_type=ParameterType.INTEGER,
                min_value=0,
                max_value=10
            )
        assert "must be >= 0" in str(excinfo.value)
        
        # 最大値超過
        with pytest.raises(ValueError) as excinfo:
            DockingParameter(
                name="above_max",
                value=11,
                parameter_type=ParameterType.INTEGER,
                min_value=0,
                max_value=10
            )
        assert "must be <= 10" in str(excinfo.value)
        
        # 境界値（最小値）
        DockingParameter(
            name="min_boundary",
            value=0,
            parameter_type=ParameterType.INTEGER,
            min_value=0,
            max_value=10
        )
        
        # 境界値（最大値）
        DockingParameter(
            name="max_boundary",
            value=10,
            parameter_type=ParameterType.INTEGER,
            min_value=0,
            max_value=10
        )
    
    def test_is_valid(self):
        """is_valid()メソッドのテスト"""
        # 有効なパラメータ
        valid_param = DockingParameter(
            name="valid_param",
            value=5,
            parameter_type=ParameterType.INTEGER,
            min_value=0,
            max_value=10
        )
        assert valid_param.is_valid() is True
        
        # 無効なパラメータはコンストラクタでエラーになるため、
        # 直接テストはできない。動的に変更する場合には注意が必要
    
    def test_to_str(self):
        """to_str()メソッドのテスト"""
        # 整数パラメータ
        int_param = DockingParameter(
            name="int_param",
            value=42,
            parameter_type=ParameterType.INTEGER
        )
        assert int_param.to_str() == "42"
        
        # 浮動小数点パラメータ
        float_param = DockingParameter(
            name="float_param",
            value=3.14,
            parameter_type=ParameterType.FLOAT
        )
        assert float_param.to_str() == "3.14"
        
        # 論理値パラメータ（小文字に変換される）
        true_param = DockingParameter(
            name="true_param",
            value=True,
            parameter_type=ParameterType.BOOLEAN
        )
        assert true_param.to_str() == "true"
        
        false_param = DockingParameter(
            name="false_param",
            value=False,
            parameter_type=ParameterType.BOOLEAN
        )
        assert false_param.to_str() == "false"
        
        # 文字列パラメータ
        str_param = DockingParameter(
            name="str_param",
            value="text",
            parameter_type=ParameterType.STRING
        )
        assert str_param.to_str() == "text"
    
    def test_to_dict(self):
        """to_dict()メソッドのテスト"""
        # 最小限のパラメータ
        min_param = DockingParameter(
            name="min_param",
            value=42,
            parameter_type=ParameterType.INTEGER
        )
        
        min_dict = min_param.to_dict()
        assert min_dict["name"] == "min_param"
        assert min_dict["value"] == 42
        assert min_dict["type"] == "INTEGER"
        assert "description" not in min_dict
        assert "min_value" not in min_dict
        assert "max_value" not in min_dict
        
        # 全オプション付きパラメータ
        full_param = DockingParameter(
            name="full_param",
            value=3.14,
            parameter_type=ParameterType.FLOAT,
            description="Pi approximation",
            min_value=0.0,
            max_value=10.0
        )
        
        full_dict = full_param.to_dict()
        assert full_dict["name"] == "full_param"
        assert full_dict["value"] == 3.14
        assert full_dict["type"] == "FLOAT"
        assert full_dict["description"] == "Pi approximation"
        assert full_dict["min_value"] == 0.0
        assert full_dict["max_value"] == 10.0


class TestDockingParameters:
    """DockingParametersクラスのテストケース"""
    
    def test_initialization_and_add(self):
        """初期化とパラメータ追加のテスト"""
        # 空のパラメータコレクション
        params = DockingParameters()
        assert len(params.parameters) == 0
        
        # パラメータの追加
        params.add(DockingParameter(
            name="exhaustiveness",
            value=8,
            parameter_type=ParameterType.INTEGER
        ))
        
        assert len(params.parameters) == 1
        assert "exhaustiveness" in params.parameters
        assert params.parameters["exhaustiveness"].value == 8
        
        # さらにパラメータを追加
        params.add(DockingParameter(
            name="energy_range",
            value=3.0,
            parameter_type=ParameterType.FLOAT
        ))
        
        assert len(params.parameters) == 2
        assert "energy_range" in params.parameters
    
    def test_get(self):
        """get()メソッドのテスト"""
        params = DockingParameters()
        
        # パラメータの追加
        params.add(DockingParameter(
            name="exhaustiveness",
            value=8,
            parameter_type=ParameterType.INTEGER
        ))
        
        # 存在するパラメータの取得
        param = params.get("exhaustiveness")
        assert param is not None
        assert param.name == "exhaustiveness"
        assert param.value == 8
        
        # 存在しないパラメータの取得
        param = params.get("non_existent")
        assert param is None
    
    def test_get_value(self):
        """get_value()メソッドのテスト"""
        params = DockingParameters()
        
        # パラメータの追加
        params.add(DockingParameter(
            name="exhaustiveness",
            value=8,
            parameter_type=ParameterType.INTEGER
        ))
        
        # 存在するパラメータの値の取得
        value = params.get_value("exhaustiveness")
        assert value == 8
        
        # 存在しないパラメータの値の取得（デフォルトなし）
        value = params.get_value("non_existent")
        assert value is None
        
        # 存在しないパラメータの値の取得（デフォルトあり）
        value = params.get_value("non_existent", default=42)
        assert value == 42
    
    def test_remove(self):
        """remove()メソッドのテスト"""
        params = DockingParameters()
        
        # パラメータの追加
        params.add(DockingParameter(
            name="param1",
            value=1,
            parameter_type=ParameterType.INTEGER
        ))
        
        params.add(DockingParameter(
            name="param2",
            value=2,
            parameter_type=ParameterType.INTEGER
        ))
        
        assert len(params.parameters) == 2
        
        # パラメータの削除
        params.remove("param1")
        
        assert len(params.parameters) == 1
        assert "param1" not in params.parameters
        assert "param2" in params.parameters
        
        # 存在しないパラメータの削除（エラーにならない）
        params.remove("non_existent")
        assert len(params.parameters) == 1
    
    def test_has(self):
        """has()メソッドのテスト"""
        params = DockingParameters()
        
        # パラメータの追加
        params.add(DockingParameter(
            name="param1",
            value=1,
            parameter_type=ParameterType.INTEGER
        ))
        
        # 存在するパラメータの確認
        assert params.has("param1") is True
        
        # 存在しないパラメータの確認
        assert params.has("non_existent") is False
    
    def test_validate(self):
        """validate()メソッドのテスト"""
        # 有効なパラメータのみを含むコレクション
        valid_params = DockingParameters()
        valid_params.add(DockingParameter(
            name="param1",
            value=1,
            parameter_type=ParameterType.INTEGER
        ))
        valid_params.add(DockingParameter(
            name="param2",
            value=True,
            parameter_type=ParameterType.BOOLEAN
        ))
        
        assert valid_params.validate() is True
        
        # 無効なパラメータを含むコレクションはコンストラクタでエラーになるため、
        # 直接テストはできない
    
    def test_to_dict(self):
        """to_dict()メソッドのテスト"""
        params = DockingParameters()
        
        # パラメータの追加
        params.add(DockingParameter(
            name="int_param",
            value=42,
            parameter_type=ParameterType.INTEGER
        ))
        
        params.add(DockingParameter(
            name="float_param",
            value=3.14,
            parameter_type=ParameterType.FLOAT
        ))
        
        params.add(DockingParameter(
            name="bool_param",
            value=True,
            parameter_type=ParameterType.BOOLEAN
        ))
        
        # 辞書への変換
        params_dict = params.to_dict()
        
        assert len(params_dict) == 3
        assert params_dict["int_param"] == 42
        assert params_dict["float_param"] == 3.14
        assert params_dict["bool_param"] is True
    
    def test_from_dict(self):
        """from_dict()メソッドのテスト"""
        # 入力辞書
        input_dict = {
            "int_param": 42,
            "float_param": 3.14,
            "bool_param": True,
            "str_param": "text"
        }
        
        # 辞書からパラメータの作成
        params = DockingParameters.from_dict(input_dict)
        
        assert len(params.parameters) == 4
        
        # 各パラメータの存在確認
        assert params.has("int_param") is True
        assert params.has("float_param") is True
        assert params.has("bool_param") is True
        assert params.has("str_param") is True
        
        # 各パラメータの値の確認
        assert params.get_value("int_param") == 42
        assert params.get_value("float_param") == 3.14
        assert params.get_value("bool_param") is True
        assert params.get_value("str_param") == "text"
        
        # 各パラメータの型の確認
        int_param = params.get("int_param")
        assert int_param is not None
        assert int_param.parameter_type == ParameterType.INTEGER
        
        float_param = params.get("float_param")
        assert float_param is not None
        assert float_param.parameter_type == ParameterType.FLOAT
        
        bool_param = params.get("bool_param")
        assert bool_param is not None
        assert bool_param.parameter_type == ParameterType.BOOLEAN
        
        str_param = params.get("str_param")
        assert str_param is not None
        assert str_param.parameter_type == ParameterType.STRING
    
    def test_merge(self):
        """merge()メソッドのテスト"""
        # 第1のパラメータコレクション
        params1 = DockingParameters()
        params1.add(DockingParameter(
            name="common_param",
            value=1,
            parameter_type=ParameterType.INTEGER
        ))
        params1.add(DockingParameter(
            name="unique_param1",
            value="unique1",
            parameter_type=ParameterType.STRING
        ))
        
        # 第2のパラメータコレクション
        params2 = DockingParameters()
        params2.add(DockingParameter(
            name="common_param",  # 同じ名前、値が異なる
            value=2,
            parameter_type=ParameterType.INTEGER
        ))
        params2.add(DockingParameter(
            name="unique_param2",
            value="unique2",
            parameter_type=ParameterType.STRING
        ))
        
        # マージ
        merged = params1.merge(params2)
        
        # 結果の確認
        assert len(merged.parameters) == 3
        
        # 共通パラメータは後のコレクション（params2）の値が優先される
        assert merged.get_value("common_param") == 2
        
        # 一意のパラメータは両方のコレクションから含まれる
        assert merged.get_value("unique_param1") == "unique1"
        assert merged.get_value("unique_param2") == "unique2"