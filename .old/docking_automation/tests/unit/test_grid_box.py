#!/usr/bin/env python3
"""
GridBoxクラスのテスト

このモジュールはドッキング計算のグリッドボックスを表すGridBoxクラスをテストします。
グリッドボックスは、ドッキング計算の探索空間を定義する重要な要素です。
"""

import pytest
from docking_automation.docking.value_object.grid_box import GridBox


class TestGridBox:
    """GridBoxクラスのテストケース"""

    def test_initialization(self):
        """初期化と基本プロパティのテスト"""
        # 通常の初期化
        grid_box = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        # 各プロパティが正しく設定されているか確認
        assert grid_box.center_x == 10.0
        assert grid_box.center_y == 20.0
        assert grid_box.center_z == 30.0
        assert grid_box.size_x == 15.0
        assert grid_box.size_y == 25.0
        assert grid_box.size_z == 35.0
    
    def test_corner_calculation(self):
        """角点（境界）計算のテスト"""
        grid_box = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=10.0,
            size_y=20.0,
            size_z=30.0
        )
        
        # 最小角点（min corner）
        min_x, min_y, min_z = grid_box.get_corner_min()
        assert min_x == 5.0  # center_x - size_x/2
        assert min_y == 10.0  # center_y - size_y/2
        assert min_z == 15.0  # center_z - size_z/2
        
        # 最大角点（max corner）
        max_x, max_y, max_z = grid_box.get_corner_max()
        assert max_x == 15.0  # center_x + size_x/2
        assert max_y == 30.0  # center_y + size_y/2
        assert max_z == 45.0  # center_z + size_z/2
    
    def test_volume_calculation(self):
        """体積計算のテスト"""
        grid_box = GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=10.0,
            size_y=10.0,
            size_z=10.0
        )
        
        # 体積計算
        volume = grid_box.get_volume()
        assert volume == 1000.0  # 10 * 10 * 10
    
    def test_validation(self):
        """バリデーションのテスト（__post_init__による検証）"""
        # 有効なグリッドボックス
        valid_grid = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        # 正常に初期化されること
        assert valid_grid is not None
        
        # 無効なグリッドボックス（サイズが負の値）
        with pytest.raises(ValueError) as excinfo:
            GridBox(
                center_x=10.0,
                center_y=20.0,
                center_z=30.0,
                size_x=-5.0,  # 負のサイズ
                size_y=25.0,
                size_z=35.0
            )
        assert "Grid box sizes must be positive" in str(excinfo.value)
        
        # 無効なグリッドボックス（サイズがゼロ）
        with pytest.raises(ValueError) as excinfo:
            GridBox(
                center_x=10.0,
                center_y=20.0,
                center_z=30.0,
                size_x=0.0,  # ゼロサイズ
                size_y=25.0,
                size_z=35.0
            )
        assert "Grid box sizes must be positive" in str(excinfo.value)
        
        # 無効なスペーシング
        with pytest.raises(ValueError) as excinfo:
            GridBox(
                center_x=10.0,
                center_y=20.0,
                center_z=30.0,
                size_x=15.0,
                size_y=25.0,
                size_z=35.0,
                spacing=-0.5  # 負のスペーシング
            )
        assert "Grid spacing must be positive" in str(excinfo.value)
    
    def test_contains_point(self):
        """点が含まれるかどうかのテスト"""
        grid_box = GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=10.0,
            size_y=10.0,
            size_z=10.0
        )
        
        # グリッドボックス内の点
        assert grid_box.contains_point(0.0, 0.0, 0.0) is True  # 中心点
        assert grid_box.contains_point(4.9, 4.9, 4.9) is True  # 境界近く（内側）
        
        # グリッドボックス外の点
        assert grid_box.contains_point(6.0, 0.0, 0.0) is False  # X軸上で外
        assert grid_box.contains_point(0.0, 6.0, 0.0) is False  # Y軸上で外
        assert grid_box.contains_point(0.0, 0.0, 6.0) is False  # Z軸上で外
        
        # 境界上の点（境界も含む場合）
        assert grid_box.contains_point(5.0, 0.0, 0.0) is True  # X軸の境界
        assert grid_box.contains_point(0.0, 5.0, 0.0) is True  # Y軸の境界
        assert grid_box.contains_point(0.0, 0.0, 5.0) is True  # Z軸の境界
    
    def test_from_dict(self):
        """辞書からの生成テスト"""
        grid_dict = {
            'center_x': 10.0,
            'center_y': 20.0,
            'center_z': 30.0,
            'size_x': 15.0,
            'size_y': 25.0,
            'size_z': 35.0
        }
        
        grid_box = GridBox.from_dict(grid_dict)
        
        assert grid_box.center_x == 10.0
        assert grid_box.center_y == 20.0
        assert grid_box.center_z == 30.0
        assert grid_box.size_x == 15.0
        assert grid_box.size_y == 25.0
        assert grid_box.size_z == 35.0
    
    def test_to_dict(self):
        """辞書への変換テスト"""
        grid_box = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        grid_dict = grid_box.to_dict()
        
        assert grid_dict['center_x'] == 10.0
        assert grid_dict['center_y'] == 20.0
        assert grid_dict['center_z'] == 30.0
        assert grid_dict['size_x'] == 15.0
        assert grid_dict['size_y'] == 25.0
        assert grid_dict['size_z'] == 35.0
        assert grid_dict['spacing'] == 1.0  # デフォルト値
    
    def test_with_padding(self):
        """パディング追加のテスト"""
        original_grid = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        # パディングを追加
        padding = 5.0
        padded_grid = original_grid.with_padding(padding)
        
        # 中心座標は変わらないことを確認
        assert padded_grid.center_x == original_grid.center_x
        assert padded_grid.center_y == original_grid.center_y
        assert padded_grid.center_z == original_grid.center_z
        
        # サイズはパディングの2倍だけ増加することを確認
        assert padded_grid.size_x == original_grid.size_x + 2 * padding
        assert padded_grid.size_y == original_grid.size_y + 2 * padding
        assert padded_grid.size_z == original_grid.size_z + 2 * padding
        
        # スペーシングは変わらないことを確認
        assert padded_grid.spacing == original_grid.spacing
        
        # 負のパディングでエラーが発生することを確認
        with pytest.raises(ValueError) as excinfo:
            original_grid.with_padding(-1.0)
        assert "Padding must be non-negative" in str(excinfo.value)
    
    def test_get_center_and_size(self):
        """中心とサイズの取得メソッドをテスト"""
        grid_box = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        # 中心を取得
        center = grid_box.get_center()
        assert center == (10.0, 20.0, 30.0)
        
        # サイズを取得
        size = grid_box.get_size()
        assert size == (15.0, 25.0, 35.0)
    
    def test_equality(self):
        """等価性のテスト"""
        grid1 = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        grid2 = GridBox(
            center_x=10.0,
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        grid3 = GridBox(
            center_x=11.0,  # 値が異なる
            center_y=20.0,
            center_z=30.0,
            size_x=15.0,
            size_y=25.0,
            size_z=35.0
        )
        
        assert grid1 == grid2  # 同じ値を持つグリッドボックスは等しい
        assert grid1 != grid3  # 異なる値を持つグリッドボックスは等しくない
        assert grid1 != "not a grid box"  # 異なる型とは等しくない