from pathlib import Path
from typing import List, Optional, Any

from openbabel import openbabel as ob
from rdkit import Chem

from docking_automation.molecule.protein import Protein
from docking_automation.molecule.compound_set import CompoundSet


class MoleculeConverter:
    """
    分子変換処理を行うクラス。
    
    OpenBabel、RDKit、meeko等を使用して、様々な分子形式間の変換を行う。
    """
    
    def protein_to_openbabel(self, protein: Protein) -> Any:
        """
        ProteinオブジェクトからOpenBabelの分子オブジェクトに変換する。
        
        Args:
            protein: 変換対象のProteinオブジェクト
            
        Returns:
            OpenBabelの分子オブジェクト
        """
        raise NotImplementedError()
    
    def openbabel_to_protein(self, ob_mol: Any, id: str, output_path: Path) -> Protein:
        """
        OpenBabelの分子オブジェクトからProteinオブジェクトに変換する。
        
        Args:
            ob_mol: 変換対象のOpenBabel分子オブジェクト
            id: 生成するProteinのID
            output_path: 出力先のパス
            
        Returns:
            生成されたProteinオブジェクト
        """
        raise NotImplementedError()
    
    def compound_to_rdkit(self, compound: CompoundSet) -> List[Any]:
        """
        CompoundSetオブジェクトからRDKitの分子オブジェクトのリストに変換する。
        
        Args:
            compound: 変換対象のCompoundSetオブジェクト
            
        Returns:
            RDKitの分子オブジェクトのリスト
        """
        raise NotImplementedError()
    
    def rdkit_to_compound(self, mols: List[Any], id: str, output_path: Path) -> CompoundSet:
        """
        RDKitの分子オブジェクトのリストからCompoundSetオブジェクトに変換する。
        
        Args:
            mols: 変換対象のRDKit分子オブジェクトのリスト
            id: 生成するCompoundSetのID
            output_path: 出力先のパス
            
        Returns:
            生成されたCompoundSetオブジェクト
        """
        raise NotImplementedError()
    
    def protein_to_pdbqt(self, protein: Protein, output_path: Path) -> Path:
        """
        Proteinオブジェクトからpdbqtファイルに変換する。
        
        Args:
            protein: 変換対象のProteinオブジェクト
            output_path: 出力先のパス
            
        Returns:
            生成されたpdbqtファイルのパス
        """
        raise NotImplementedError()
    
    def compound_to_pdbqt(self, compound: CompoundSet, output_path: Path) -> Path:
        """
        CompoundSetオブジェクトからpdbqtファイルに変換する。
        
        Args:
            compound: 変換対象のCompoundSetオブジェクト
            output_path: 出力先のパス
            
        Returns:
            生成されたpdbqtファイルのパス
        """
        raise NotImplementedError()
    
    def pdbqt_to_sdf(self, pdbqt_path: Path, output_path: Path) -> Path:
        """
        pdbqtファイルからsdfファイルに変換する。
        
        Args:
            pdbqt_path: 変換対象のpdbqtファイルのパス
            output_path: 出力先のパス
            
        Returns:
            生成されたsdfファイルのパス
        """
        raise NotImplementedError()
    
    def pdbqt_to_rdkit(self, pdbqt_path: Path) -> List[Any]:
        """
        pdbqtファイルからRDKitの分子オブジェクトのリストに変換する。
        
        Args:
            pdbqt_path: 変換対象のpdbqtファイルのパス
            
        Returns:
            RDKitの分子オブジェクトのリスト
        """
        raise NotImplementedError()