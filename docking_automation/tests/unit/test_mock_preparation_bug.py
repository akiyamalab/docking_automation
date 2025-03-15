import tempfile
import os
import pytest
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.infrastructure.tools.mock.mock_molecule_preparation_service import MockMoleculePreparationService

def test_compound_has_no_is_prepared_attribute():
    """Compoundクラスにis_prepared属性が存在しないことを確認するテスト"""
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp.write(b"ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00\n")
        tmp_path = tmp.name
    
    try:
        # Compoundインスタンスを作成
        compound = Compound(
            id="test",
            path=tmp_path,
            structure=None,
            format=MoleculeFormat(type=FormatType.PDB),
            metadata={}
        )
        
        # is_prepared属性が存在しないことを確認
        assert not hasattr(compound, 'is_prepared')
        
        # format.typeアクセスは可能
        assert compound.format is not None
        assert compound.format.type == FormatType.PDB
        
        # propertiesは存在しない
        assert not hasattr(compound, 'properties')
        
        # preparation_methodは存在しない
        assert not hasattr(compound, 'preparation_method')
    finally:
        # 一時ファイルを削除
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def test_protein_has_no_special_attributes():
    """Proteinクラスに特殊な属性が存在しないことを確認するテスト"""
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp.write(b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\n")
        tmp_path = tmp.name
    
    try:
        # Proteinインスタンスを作成
        protein = Protein(
            id="test",
            path=tmp_path,
            structure=None,
            format=MoleculeFormat(type=FormatType.PDB),
            chains={"A"},
            metadata={}
        )
        
        # 特殊な属性が存在しないことを確認
        assert not hasattr(protein, 'is_prepared')
        assert not hasattr(protein, 'has_water')
        assert not hasattr(protein, 'has_hydrogens')
        assert not hasattr(protein, 'properties')
        assert not hasattr(protein, 'preparation_method')
        assert not hasattr(protein, 'active_site_residues')
        
        # format.typeアクセスは可能
        assert protein.format is not None
        assert protein.format.type == FormatType.PDB
        
        # chainsは存在する
        assert hasattr(protein, 'chains')
        assert "A" in protein.chains
    finally:
        # 一時ファイルを削除
        if os.path.exists(tmp_path):
            os.unlink(tmp_path)


def test_mock_preparation_service_bug():
    """MockMoleculePreparationServiceの問題を確認するテスト - バグ修正後の正常動作を確認"""
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp.write(b"ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00\n")
        compound_path = tmp.name
        
    with tempfile.NamedTemporaryFile(suffix='.pdb', delete=False) as tmp:
        tmp.write(b"ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\n")
        protein_path = tmp.name
    
    try:
        # テスト対象のインスタンスを作成
        compound = Compound(
            id="test_compound",
            path=compound_path,
            structure=None,
            format=MoleculeFormat(type=FormatType.PDB),
            metadata={}
        )
        
        protein = Protein(
            id="test_protein",
            path=protein_path,
            structure=None,
            format=MoleculeFormat(type=FormatType.PDB),
            chains={"A"},
            metadata={}
        )
        
        # MockMoleculePreparationServiceを作成
        service = MockMoleculePreparationService()
        
        # 修正後のコードでは、is_prepared属性がなくてもエラーが発生しないことを確認
        try:
            # 準備メソッドを呼び出す
            prepared_compound = service.prepare_ligand(compound)
            prepared_protein = service.prepare_receptor(protein)
            
            # 返されたオブジェクトが正しい型であることを確認
            assert isinstance(prepared_compound, Compound)
            assert isinstance(prepared_protein, Protein)
            
            # 返されたオブジェクトにメタデータが設定されていることを確認
            assert "is_prepared" in prepared_compound.metadata
            assert prepared_compound.metadata["is_prepared"] is True
            
            assert "is_prepared" in prepared_protein.metadata
            assert prepared_protein.metadata["is_prepared"] is True
            
        except AttributeError as e:
            pytest.fail(f"AttributeError was raised: {e}")
    finally:
        # 一時ファイルを削除
        if os.path.exists(compound_path):
            os.unlink(compound_path)
        if os.path.exists(protein_path):
            os.unlink(protein_path)