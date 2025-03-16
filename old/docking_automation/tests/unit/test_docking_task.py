#!/usr/bin/env python3
"""
DockingTaskクラスのテスト

このモジュールはドッキング計算タスクを表すDockingTaskクラスをテストします。
DockingTaskはリガンド、レセプター、設定を含み、ドッキング計算のコンテキストを表します。
"""

import time
import pytest
from unittest.mock import patch, MagicMock
from dataclasses import dataclass, field
from typing import Dict, Any, Optional, List, Set, Tuple

from docking_automation.docking.entity.docking_task import DockingTask, TaskStatus
from docking_automation.docking.entity.ligand import Ligand
from docking_automation.docking.entity.receptor import Receptor
from docking_automation.docking.value_object.docking_configuration import DockingConfiguration
from docking_automation.docking.value_object.grid_box import GridBox
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein
from docking_automation.docking.value_object.docking_parameter import DockingParameters


class TestDockingTask:
    """DockingTaskクラスのテストケース"""
    
    @pytest.fixture
    def mock_compound(self):
        """テスト用のモック化合物"""
        compound = MagicMock(spec=Compound)
        compound.id = "C001"
        compound.name = "test_compound"
        compound.path = "test/compound.sdf"
        compound.format.type = MagicMock()  # Ligandのis_preparedメソッドでチェックされる
        compound.is_prepared = True  # Ligandのis_preparedメソッドでチェックされる
        return compound
    
    @pytest.fixture
    def mock_protein(self):
        """テスト用のモックタンパク質"""
        protein = MagicMock(spec=Protein)
        protein.id = "P001"
        protein.name = "test_protein"
        protein.path = "test/protein.pdb"
        protein.chains = {"A", "B"}
        return protein
    
    @pytest.fixture
    def mock_ligand(self, mock_compound):
        """テスト用のモックリガンド"""
        ligand = MagicMock(spec=Ligand)
        ligand.compound = mock_compound
        ligand.name = mock_compound.name
        ligand.is_prepared.return_value = True
        return ligand
    
    @pytest.fixture
    def mock_receptor(self, mock_protein):
        """テスト用のモックレセプター"""
        receptor = MagicMock(spec=Receptor)
        receptor.protein = mock_protein
        receptor.name = mock_protein.name
        receptor.is_prepared.return_value = True
        return receptor
    
    @pytest.fixture
    def mock_grid_box(self):
        """テスト用のモックグリッドボックス"""
        return GridBox(
            center_x=0.0,
            center_y=0.0,
            center_z=0.0,
            size_x=20.0,
            size_y=20.0,
            size_z=20.0
        )
    
    @pytest.fixture
    def mock_configuration(self, mock_grid_box):
        """テスト用のモック設定"""
        config = MagicMock(spec=DockingConfiguration)
        config.grid_box = mock_grid_box
        config.parameters = DockingParameters()
        config.validate.return_value = True
        return config
    
    def test_initialization(self, mock_ligand, mock_receptor, mock_configuration):
        """初期化と基本プロパティのテスト"""
        # 明示的なIDを指定
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        assert task.id == "T001"
        assert task.ligand is mock_ligand
        assert task.receptor is mock_receptor
        assert task.configuration is mock_configuration
        assert task.status == TaskStatus.PENDING
        assert task.error_message is None
        assert task.metadata == {}
        assert task.created_at is not None
        assert task.updated_at is not None
        assert task.timeout_seconds == 3600  # デフォルト値
        assert task.molecular_weight is None
        
        # ID自動生成テスト
        task = DockingTask(
            id="",  # 空文字列の場合、UUIDが自動生成される
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        assert task.id != ""  # UUIDが生成されている
        assert len(task.id) == 32  # UUIDのhex文字列は32文字
    
    def test_update_status(self, mock_ligand, mock_receptor, mock_configuration):
        """ステータス更新のテスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # 初期状態
        assert task.status == TaskStatus.PENDING
        assert task.error_message is None
        
        # ステータスを更新
        original_updated_at = task.updated_at
        assert original_updated_at is not None  # Noneでないことを確認
        time.sleep(0.001)  # 微小な時間待機
        
        task.update_status(TaskStatus.PREPARING)
        assert task.status == TaskStatus.PREPARING
        assert task.error_message is None
        assert task.updated_at is not None
        assert task.updated_at > original_updated_at  # 更新時刻が更新されている
        
        # エラーメッセージ付きでステータスを更新
        original_updated_at = task.updated_at
        assert original_updated_at is not None  # Noneでないことを確認
        time.sleep(0.001)  # 微小な時間待機
        
        task.update_status(TaskStatus.FAILED, "Something went wrong")
        assert task.status == TaskStatus.FAILED
        assert task.error_message == "Something went wrong"
        assert task.updated_at is not None
        assert task.updated_at > original_updated_at  # 更新時刻が更新されている
    
    def test_get_receptors_single(self, mock_ligand, mock_receptor, mock_configuration):
        """単一レセプターの取得テスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,  # 単一レセプター
            configuration=mock_configuration
        )
        
        receptors = task.get_receptors()
        assert isinstance(receptors, list)
        assert len(receptors) == 1
        assert receptors[0] is mock_receptor
    
    def test_get_receptors_multiple(self, mock_ligand, mock_receptor, mock_configuration):
        """複数レセプターの取得テスト"""
        # 追加のモックレセプター
        mock_receptor2 = MagicMock(spec=Receptor)
        mock_receptor2.name = "test_receptor2"
        mock_receptor2.is_prepared.return_value = True
        
        # 複数レセプターでタスクを作成
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=[mock_receptor, mock_receptor2],  # レセプターのリスト
            configuration=mock_configuration
        )
        
        receptors = task.get_receptors()
        assert isinstance(receptors, list)
        assert len(receptors) == 2
        assert receptors[0] is mock_receptor
        assert receptors[1] is mock_receptor2
    
    def test_is_multi_receptor(self, mock_ligand, mock_receptor, mock_configuration):
        """複数レセプター判定のテスト"""
        # 単一レセプター
        task_single = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        assert task_single.is_multi_receptor() is False
        
        # 複数レセプター
        mock_receptor2 = MagicMock(spec=Receptor)
        
        task_multi = DockingTask(
            id="T002",
            ligand=mock_ligand,
            receptor=[mock_receptor, mock_receptor2],
            configuration=mock_configuration
        )
        
        assert task_multi.is_multi_receptor() is True
        
        # 単一要素のリスト（マルチレセプターではない）
        task_single_list = DockingTask(
            id="T003",
            ligand=mock_ligand,
            receptor=[mock_receptor],
            configuration=mock_configuration
        )
        
        assert task_single_list.is_multi_receptor() is False
    
    def test_get_primary_receptor(self, mock_ligand, mock_receptor, mock_configuration):
        """主要レセプター取得のテスト"""
        # 単一レセプター
        task_single = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        primary = task_single.get_primary_receptor()
        assert primary is mock_receptor
        
        # 複数レセプター
        mock_receptor2 = MagicMock(spec=Receptor)
        
        task_multi = DockingTask(
            id="T002",
            ligand=mock_ligand,
            receptor=[mock_receptor, mock_receptor2],
            configuration=mock_configuration
        )
        
        primary = task_multi.get_primary_receptor()
        assert primary is mock_receptor  # 最初のレセプターが返される
        
        # 空のレセプターリスト
        task_empty = DockingTask(
            id="T003",
            ligand=mock_ligand,
            receptor=[],
            configuration=mock_configuration
        )
        
        with pytest.raises(ValueError) as excinfo:
            task_empty.get_primary_receptor()
        assert "No receptors available" in str(excinfo.value)
    
    def test_is_ready(self, mock_ligand, mock_receptor, mock_configuration):
        """実行準備完了確認のテスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # 準備完了の場合
        mock_ligand.is_prepared.return_value = True
        mock_receptor.is_prepared.return_value = True
        mock_configuration.validate.return_value = True
        
        assert task.is_ready() is True
        
        # リガンドが準備されていない場合
        mock_ligand.is_prepared.return_value = False
        mock_receptor.is_prepared.return_value = True
        mock_configuration.validate.return_value = True
        
        assert task.is_ready() is False
        
        # レセプターが準備されていない場合
        mock_ligand.is_prepared.return_value = True
        mock_receptor.is_prepared.return_value = False
        mock_configuration.validate.return_value = True
        
        assert task.is_ready() is False
        
        # 設定が無効な場合
        mock_ligand.is_prepared.return_value = True
        mock_receptor.is_prepared.return_value = True
        mock_configuration.validate.return_value = False
        
        assert task.is_ready() is False
        
        # 分子量フィルターテスト
        mock_ligand.is_prepared.return_value = True
        mock_receptor.is_prepared.return_value = True
        mock_configuration.validate.return_value = True
        
        # 分子量が制限未満
        task.molecular_weight = 699
        assert task.is_ready() is True
        
        # 分子量が制限以上
        task.molecular_weight = 700
        assert task.is_ready() is False
    
    def test_prepare_for_execution(self, mock_ligand, mock_receptor, mock_configuration):
        """実行準備のテスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # 準備完了の場合
        mock_ligand.is_prepared.return_value = True
        mock_receptor.is_prepared.return_value = True
        mock_configuration.validate.return_value = True
        
        result = task.prepare_for_execution()
        assert result is True
        assert task.status == TaskStatus.READY
        
        # 既にREADYの場合
        task.status = TaskStatus.READY
        result = task.prepare_for_execution()
        assert result is True  # 何もせずにTrueを返す
        
        # 既にRUNNINGの場合
        task.status = TaskStatus.RUNNING
        result = task.prepare_for_execution()
        assert result is True  # 何もせずにTrueを返す
        
        # 既にCOMPLETEDの場合
        task.status = TaskStatus.COMPLETED
        result = task.prepare_for_execution()
        assert result is True  # 何もせずにTrueを返す
        
        # 準備できない場合
        task.status = TaskStatus.PENDING
        mock_ligand.is_prepared.return_value = False  # リガンドが準備されていない
        
        result = task.prepare_for_execution()
        assert result is False
        assert task.status == TaskStatus.FAILED
        assert task.error_message is not None
        assert "not properly prepared" in task.error_message
    
    def test_status_marking_methods(self, mock_ligand, mock_receptor, mock_configuration):
        """ステータスマーキングメソッドのテスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # 実行中にマーク
        task.mark_as_running()
        assert task.status == TaskStatus.RUNNING
        
        # 完了にマーク
        task.mark_as_completed()
        assert task.status == TaskStatus.COMPLETED
        
        # 失敗にマーク
        task.mark_as_failed("Test error")
        assert task.status == TaskStatus.FAILED
        assert task.error_message == "Test error"
        
        # キャンセルにマーク
        task.mark_as_cancelled()
        assert task.status == TaskStatus.CANCELLED
    
    def test_status_check_methods(self, mock_ligand, mock_receptor, mock_configuration):
        """ステータスチェックメソッドのテスト"""
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # PENDINGステータス
        task.status = TaskStatus.PENDING
        assert task.is_active() is False
        assert task.is_finished() is False
        
        # PREPARINGステータス
        task.status = TaskStatus.PREPARING
        assert task.is_active() is True
        assert task.is_finished() is False
        
        # READYステータス
        task.status = TaskStatus.READY
        assert task.is_active() is True
        assert task.is_finished() is False
        
        # RUNNINGステータス
        task.status = TaskStatus.RUNNING
        assert task.is_active() is True
        assert task.is_finished() is False
        
        # COMPLETEDステータス
        task.status = TaskStatus.COMPLETED
        assert task.is_active() is False
        assert task.is_finished() is True
        
        # FAILEDステータス
        task.status = TaskStatus.FAILED
        assert task.is_active() is False
        assert task.is_finished() is True
        
        # CANCELLEDステータス
        task.status = TaskStatus.CANCELLED
        assert task.is_active() is False
        assert task.is_finished() is True
    
    def test_get_time_elapsed(self, mock_ligand, mock_receptor, mock_configuration):
        """経過時間取得のテスト"""
        # 現在時刻をモックする
        with patch('time.time') as mock_time:
            # 初期時刻を設定
            start_time = 1000.0
            mock_time.return_value = start_time
            
            task = DockingTask(
                id="T001",
                ligand=mock_ligand,
                receptor=mock_receptor,
                configuration=mock_configuration
            )
            
            # created_atが設定されていることを確認
            assert task.created_at == start_time
            
            # 進行中のタスクの経過時間
            mock_time.return_value = start_time + 10.0  # 10秒後
            elapsed = task.get_time_elapsed()
            assert elapsed == 10.0
            
            # 完了したタスクの経過時間
            task.mark_as_completed()  # updated_atが更新される
            mock_time.return_value = start_time + 20.0  # さらに10秒後
            
            # updated_atの時点での時間が使われる
            assert task.get_time_elapsed() == 10.0
            
            # created_atがNoneの場合
            task.created_at = None
            assert task.get_time_elapsed() is None
    
    def test_is_timed_out(self, mock_ligand, mock_receptor, mock_configuration):
        """タイムアウト確認のテスト"""
        # 現在時刻をモックする
        with patch('time.time') as mock_time:
            # 初期時刻を設定
            start_time = 1000.0
            mock_time.return_value = start_time
            
            task = DockingTask(
                id="T001",
                ligand=mock_ligand,
                receptor=mock_receptor,
                configuration=mock_configuration,
                timeout_seconds=30  # 30秒タイムアウト
            )
            
            # タイムアウトしていない場合
            mock_time.return_value = start_time + 20.0  # 20秒後
            assert task.is_timed_out() is False
            
            # タイムアウトした場合
            mock_time.return_value = start_time + 31.0  # 31秒後
            assert task.is_timed_out() is True
            
            # created_atがNoneの場合
            task.created_at = None
            assert task.is_timed_out() is False
    
    def test_str_representation_single_receptor(self, mock_ligand, mock_receptor, mock_configuration):
        """単一レセプターの文字列表現テスト"""
        mock_receptor.name = "TestReceptor"
        mock_ligand.name = "TestLigand"
        
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=mock_receptor,
            configuration=mock_configuration
        )
        
        # 文字列表現
        str_rep = str(task)
        assert "T001" in str_rep
        assert "TestLigand" in str_rep
        assert "TestReceptor" in str_rep
        assert "PENDING" in str_rep
    
    def test_str_representation_multiple_receptors(self, mock_ligand, mock_configuration):
        """複数レセプターの文字列表現テスト"""
        mock_receptor1 = MagicMock(spec=Receptor)
        mock_receptor1.name = "TestReceptor1"
        
        mock_receptor2 = MagicMock(spec=Receptor)
        mock_receptor2.name = "TestReceptor2"
        
        mock_receptor3 = MagicMock(spec=Receptor)
        mock_receptor3.name = "TestReceptor3"
        
        mock_ligand.name = "TestLigand"
        
        task = DockingTask(
            id="T001",
            ligand=mock_ligand,
            receptor=[mock_receptor1, mock_receptor2, mock_receptor3],
            configuration=mock_configuration
        )
        
        # 文字列表現
        str_rep = str(task)
        assert "T001" in str_rep
        assert "TestLigand" in str_rep
        assert "TestReceptor1" in str_rep  # 主要レセプター名
        assert "+2 more" in str_rep  # 追加レセプター数
        assert "PENDING" in str_rep