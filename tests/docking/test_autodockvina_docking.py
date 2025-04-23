import os
import tempfile
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from docking_automation.docking import DockingResult
from docking_automation.docking.autodockvina_docking import AutoDockVina, AutoDockVinaParameters
from docking_automation.docking.docking_parameters import CommonDockingParameters, DockingParameters
from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.compound_set import CompoundSet
from docking_automation.molecule.protein import Protein


class TestDockingResult:
    """DockingResultクラスのテスト"""

    def test_docking_result_has_hash_attributes(self):
        """DockingResultオブジェクトがhash属性を持つことを確認する"""
        # テスト用のパラメータ
        result_path = Path("/tmp/test.sdf")
        protein_id = "test_protein"
        compound_set_id = "test_compound_set"
        compound_index = 0
        docking_score = -8.0
        protein_content_hash = "test_protein_hash"
        compoundset_content_hash = "test_compound_hash"
        metadata = {"test": "metadata"}

        # DockingResultオブジェクトを作成
        result = DockingResult(
            result_path=result_path,
            protein_id=protein_id,
            compound_set_id=compound_set_id,
            compound_index=compound_index,
            docking_score=docking_score,
            protein_content_hash=protein_content_hash,
            compoundset_content_hash=compoundset_content_hash,
            metadata=metadata,
        )

        # 検証
        assert hasattr(result, "protein_content_hash")
        assert hasattr(result, "compoundset_content_hash")
        assert result.protein_content_hash == protein_content_hash
        assert result.compoundset_content_hash == compoundset_content_hash


@pytest.fixture
def mock_vina():
    """Vinaクラスのモックを作成する"""
    with patch("docking_automation.docking.autodockvina_docking.Vina") as mock_vina:
        # モックVinaインスタンスの設定
        mock_vina_instance = MagicMock()
        mock_vina_instance.energies.return_value = [[[-8.0]]]
        mock_vina.return_value = mock_vina_instance
        yield mock_vina


@pytest.fixture
def mock_converter():
    """MoleculeConverterクラスのモックを作成する"""
    with patch("docking_automation.docking.autodockvina_docking.MoleculeConverter") as mock_converter:
        # モックConverterインスタンスの設定
        mock_converter_instance = MagicMock()
        mock_converter.return_value = mock_converter_instance
        yield mock_converter


@pytest.fixture
def mock_protein():
    """Proteinクラスのモックを作成する"""
    mock_protein = MagicMock(spec=Protein)
    mock_protein.id = "test_protein"
    mock_protein.content_hash = "test_protein_hash"
    
    # PreprocessedProteinのモック
    mock_preprocessed = MagicMock()
    mock_preprocessed.file_path = Path("/tmp/test_protein.pdbqt")
    mock_preprocessed.content_hash = "test_protein_hash"
    
    return mock_preprocessed


@pytest.fixture
def mock_compound_set():
    """CompoundSetクラスのモックを作成する"""
    mock_compound_set = MagicMock(spec=CompoundSet)
    mock_compound_set.id = "test_compound_set"
    mock_compound_set.content_hash = "test_compound_hash"
    
    # PreprocessedCompoundSetのモック
    mock_preprocessed = MagicMock()
    mock_preprocessed.file_paths = [Path("/tmp/compound_0.pdbqt"), Path("/tmp/compound_1.pdbqt")]
    mock_preprocessed.content_hash = "test_compound_hash"
    mock_preprocessed.get_compound_hash.return_value = "test_compound_hash"
    
    # プロパティの設定
    properties = {
        "id": "test_compound_set",
        "path": "/tmp/test_compounds.sdf",
        "compound_count": 2,
        "content_hash": "test_compound_hash",
    }
    mock_preprocessed.get_properties.return_value = properties
    
    return mock_preprocessed


@pytest.fixture
def mock_grid_box():
    """GridBoxクラスのモックを作成する"""
    mock_grid_box = MagicMock(spec=GridBox)
    mock_grid_box.center = [0.0, 0.0, 0.0]
    mock_grid_box.size = [20.0, 20.0, 20.0]
    return mock_grid_box


class TestAutoDockVina:
    """AutoDockVinaクラスのテスト"""

    def test_dock_with_indices(self, mock_vina, mock_converter, mock_protein, mock_compound_set, mock_grid_box, tmp_path):
        """インデックスリストが設定されている場合のdockメソッドのテスト"""
        # プロパティにindicesを追加
        properties = mock_compound_set.get_properties()
        properties["indices"] = [0, 2]
        mock_compound_set.get_properties.return_value = properties
        
        # AutoDockVinaインスタンスを作成
        vina = AutoDockVina()
        vina.converter = mock_converter
        
        # パラメータを設定
        common_params = CommonDockingParameters(
            protein=mock_protein, compound_set=mock_compound_set, grid_box=mock_grid_box
        )
        specific_params = AutoDockVinaParameters(exhaustiveness=1, num_modes=1)
        parameters = DockingParameters(common=common_params, specific=specific_params)
        
        # dockメソッドを実行
        results = vina.dock(parameters)
        
        # 結果を検証
        assert len(results) == 2
        assert results[0].compound_index == 0
        assert results[1].compound_index == 2

    def test_dock_with_index_range(self, mock_vina, mock_converter, mock_protein, mock_compound_set, mock_grid_box, tmp_path):
        """インデックス範囲が設定されている場合のdockメソッドのテスト"""
        # プロパティにindex_rangeを追加
        properties = mock_compound_set.get_properties()
        properties["index_range"] = {"start": 10, "end": 12, "total_compounds": 20}
        mock_compound_set.get_properties.return_value = properties
        
        # AutoDockVinaインスタンスを作成
        vina = AutoDockVina()
        vina.converter = mock_converter
        
        # パラメータを設定
        common_params = CommonDockingParameters(
            protein=mock_protein, compound_set=mock_compound_set, grid_box=mock_grid_box
        )
        specific_params = AutoDockVinaParameters(exhaustiveness=1, num_modes=1)
        parameters = DockingParameters(common=common_params, specific=specific_params)
        
        # dockメソッドを実行
        results = vina.dock(parameters)
        
        # 結果を検証
        assert len(results) == 2
        assert results[0].compound_index == 10
        assert results[1].compound_index == 11

    def test_dock_without_indices(self, mock_vina, mock_converter, mock_protein, mock_compound_set, mock_grid_box, tmp_path):
        """インデックスが設定されていない場合のdockメソッドのテスト"""
        # AutoDockVinaインスタンスを作成
        vina = AutoDockVina()
        vina.converter = mock_converter
        
        # パラメータを設定
        common_params = CommonDockingParameters(
            protein=mock_protein, compound_set=mock_compound_set, grid_box=mock_grid_box
        )
        specific_params = AutoDockVinaParameters(exhaustiveness=1, num_modes=1)
        parameters = DockingParameters(common=common_params, specific=specific_params)
        
        # dockメソッドを実行
        results = vina.dock(parameters)
        
        # 結果を検証
        assert len(results) == 2
        assert results[0].compound_index == 0
        assert results[1].compound_index == 1

    def test_run_docking_with_compound_indices(self, mock_vina, mock_converter, mock_protein, mock_compound_set, mock_grid_box, tmp_path):
        """compound_indicesパラメータを指定した場合のrun_dockingメソッドのテスト"""
        # AutoDockVinaインスタンスを作成
        vina = AutoDockVina()
        vina.converter = mock_converter
        
        # with_indicesメソッドのモックを作成
        original_compound_set = MagicMock(spec=CompoundSet)
        original_compound_set.with_indices.return_value = mock_compound_set
        
        # _preprocess_proteinと_preprocess_compound_setメソッドをモック
        vina._preprocess_protein = MagicMock(return_value=mock_protein)
        vina._preprocess_compound_set = MagicMock(return_value=mock_compound_set)
        
        # dockメソッドをモック
        vina.dock = MagicMock(return_value=[
            MagicMock(spec=DockingResult),
            MagicMock(spec=DockingResult)
        ])
        
        # run_dockingメソッドを実行
        result_collection = vina.run_docking(
            protein=MagicMock(spec=Protein),
            compound_set=original_compound_set,
            grid_box=mock_grid_box,
            additional_params=AutoDockVinaParameters(),
            compound_indices={0, 2}
        )
        
        # 結果を検証
        assert original_compound_set.with_indices.called
        assert original_compound_set.with_indices.call_args[0][0] == {0, 2}
        assert len(result_collection) == 2

    def test_run_docking_without_compound_indices(self, mock_vina, mock_converter, mock_protein, mock_compound_set, mock_grid_box, tmp_path):
        """compound_indicesパラメータを指定しない場合のrun_dockingメソッドのテスト"""
        # AutoDockVinaインスタンスを作成
        vina = AutoDockVina()
        vina.converter = mock_converter
        
        # _preprocess_proteinと_preprocess_compound_setメソッドをモック
        vina._preprocess_protein = MagicMock(return_value=mock_protein)
        vina._preprocess_compound_set = MagicMock(return_value=mock_compound_set)
        
        # dockメソッドをモック
        vina.dock = MagicMock(return_value=[
            MagicMock(spec=DockingResult),
            MagicMock(spec=DockingResult)
        ])
        
        # run_dockingメソッドを実行
        result_collection = vina.run_docking(
            protein=MagicMock(spec=Protein),
            compound_set=mock_compound_set,
            grid_box=mock_grid_box,
            additional_params=AutoDockVinaParameters()
        )
        
        # 結果を検証
        assert len(result_collection) == 2
