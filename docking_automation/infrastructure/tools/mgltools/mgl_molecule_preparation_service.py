import os
import tempfile
import subprocess
import logging
from typing import List, Optional, Dict, Any, Union, cast, TypeVar

from docking_automation.domain.molecule.service.molecule_preparation_service import MoleculePreparationService
from docking_automation.domain.molecule.entity.compound import Compound
from docking_automation.domain.molecule.entity.protein import Protein

# 分子型のためのUnion型
MoleculeUnion = Union[Compound, Protein]
from docking_automation.domain.molecule.value_object.molecule_property import MoleculeProperty
from docking_automation.domain.molecule.value_object.molecule_format import FormatType, MoleculeFormat

# ロガーの設定
logger = logging.getLogger(__name__)


class MGLMoleculePreparationService(MoleculePreparationService):
    """MGLToolsを使用した分子準備サービスの実装"""
    
    def __init__(
        self,
        prepare_ligand_path: str = "prepare_ligand",
        prepare_receptor_path: str = "prepare_receptor",
        obabel_path: str = "obabel"
    ):
        """コンストラクタ
        
        Args:
            prepare_ligand_path: prepare_ligandコマンドの実行パス
            prepare_receptor_path: prepare_receptorコマンドの実行パス
            obabel_path: OpenBabelの実行ファイルパス（形式変換に使用）
        """
        self.prepare_ligand_path = prepare_ligand_path
        self.prepare_receptor_path = prepare_receptor_path
        self.obabel_path = obabel_path
    
    def prepare_ligand(self, compound: Compound, method: Optional[str] = None) -> Compound:
        """化合物をリガンドとして準備する
        
        典型的な準備手順：
        1. 3D構造の生成／最適化
        2. 水素原子の追加
        3. 電荷の計算と割り当て
        4. PDBQTフォーマットへの変換
        
        Args:
            compound: 準備する化合物
            method: 準備方法（"gasteiger", "am1bcc"など、実装によって異なる）
            
        Returns:
            準備された化合物
            
        Raises:
            ValueError: 準備中に無効なパラメータが検出された場合
            RuntimeError: 準備処理に失敗した場合
        """
        # 既に準備済みの場合はそのまま返す
        if compound.is_prepared and compound.format.type == FormatType.PDBQT:
            return compound
        
        # 準備するメソッドの決定
        prep_method = method or "gasteiger"  # デフォルトはGasteiger電荷
        
        try:
            # 一時ディレクトリの作成
            with tempfile.TemporaryDirectory() as temp_dir:
                # 入力ファイルの準備
                input_file = compound.path
                input_created = False
                
                # 入力ファイルがない場合は一時ファイルを作成
                if not input_file or not os.path.exists(input_file):
                    # 一時的なPDBファイルを作成
                    input_file = os.path.join(temp_dir, f"{compound.id}.pdb")
                    input_created = True
                    # ファイルに構造を書き込む処理（プレースホルダー）
                    with open(input_file, 'w') as f:
                        f.write("ATOM      1  C   LIG A   1       0.000   0.000   0.000  1.00  0.00\n")
                
                # PDB形式でない場合はOpenBabelで変換
                if compound.format.type != FormatType.PDB:
                    # OpenBabelを使ってPDBに変換
                    pdb_file = os.path.join(temp_dir, f"{compound.id}.pdb")
                    if not self._convert_to_pdb(input_file, pdb_file):
                        raise RuntimeError(f"Failed to convert {compound.id} to PDB format")
                    input_file = pdb_file
                
                # 出力ファイルパスの設定
                output_file = os.path.join(temp_dir, f"{compound.id}.pdbqt")
                
                # prepare_ligandコマンドの構築
                cmd = [
                    self.prepare_ligand_path,
                    "-l", input_file,
                    "-o", output_file,
                    "-A", "hydrogens"  # 水素を追加
                ]
                
                # 電荷計算方法に応じたオプション
                if prep_method == "gasteiger":
                    cmd.extend(["-C"])  # Gasteiger電荷
                
                # コマンドを実行
                logger.info(f"Running prepare_ligand command: {' '.join(cmd)}")
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                # 結果を確認
                if result.returncode != 0:
                    logger.error(f"prepare_ligand failed: {result.stderr}")
                    raise RuntimeError(f"prepare_ligand failed: {result.stderr}")
                
                # 出力ファイルが存在するか確認
                if not os.path.exists(output_file):
                    raise RuntimeError("Output file was not created")
                
                # 新しい化合物オブジェクトを作成
                # 実際には、新しい構造データを読み込む必要があります
                prepared_compound = Compound(
                    id=compound.id,
                    structure=compound.structure,  # 本来は変換後の構造を読み込む
                    format=MoleculeFormat(type=FormatType.PDBQT),
                    properties=compound.properties.copy(),
                    path=output_file,  # 一時ファイルのパス（後で永続化が必要）
                    metadata=compound.metadata.copy(),
                    is_prepared=True,
                    preparation_method=prep_method
                )
                
                return prepared_compound
                
        except Exception as e:
            logger.error(f"Error preparing ligand: {e}")
            raise RuntimeError(f"Error preparing ligand: {str(e)}")
    
    def prepare_receptor(self, protein: Protein, method: Optional[str] = None) -> Protein:
        """タンパク質をレセプターとして準備する
        
        典型的な準備手順：
        1. 水分子の削除（オプション）
        2. 水素原子の追加
        3. 電荷の計算と割り当て
        4. PDBQTフォーマットへの変換
        
        Args:
            protein: 準備するタンパク質
            method: 準備方法（実装によって異なる）
            
        Returns:
            準備されたタンパク質
            
        Raises:
            ValueError: 準備中に無効なパラメータが検出された場合
            RuntimeError: 準備処理に失敗した場合
        """
        # 既に準備済みの場合はそのまま返す
        if protein.is_prepared and protein.format.type == FormatType.PDBQT:
            return protein
        
        # 準備するメソッドの決定
        prep_method = method or "default"
        
        try:
            # 一時ディレクトリの作成
            with tempfile.TemporaryDirectory() as temp_dir:
                # 入力ファイルの準備
                input_file = protein.path
                input_created = False
                
                # 入力ファイルがない場合は一時ファイルを作成
                if not input_file or not os.path.exists(input_file):
                    # 一時的なPDBファイルを作成
                    input_file = os.path.join(temp_dir, f"{protein.id}.pdb")
                    input_created = True
                    # ファイルに構造を書き込む処理（プレースホルダー）
                    with open(input_file, 'w') as f:
                        f.write("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00  0.00\n")
                
                # PDB形式でない場合はOpenBabelで変換
                if protein.format.type != FormatType.PDB:
                    # OpenBabelを使ってPDBに変換
                    pdb_file = os.path.join(temp_dir, f"{protein.id}.pdb")
                    if not self._convert_to_pdb(input_file, pdb_file):
                        raise RuntimeError(f"Failed to convert {protein.id} to PDB format")
                    input_file = pdb_file
                
                # 水分子を削除する場合
                if not protein.has_water:
                    no_water_file = os.path.join(temp_dir, f"{protein.id}_nowater.pdb")
                    if not self._remove_water(input_file, no_water_file):
                        raise RuntimeError("Failed to remove water molecules")
                    input_file = no_water_file
                
                # 出力ファイルパスの設定
                output_file = os.path.join(temp_dir, f"{protein.id}.pdbqt")
                
                # prepare_receptorコマンドの構築
                cmd = [
                    self.prepare_receptor_path,
                    "-r", input_file,
                    "-o", output_file,
                    "-A", "hydrogens"  # 水素を追加
                ]
                
                # コマンドを実行
                logger.info(f"Running prepare_receptor command: {' '.join(cmd)}")
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                # 結果を確認
                if result.returncode != 0:
                    logger.error(f"prepare_receptor failed: {result.stderr}")
                    raise RuntimeError(f"prepare_receptor failed: {result.stderr}")
                
                # 出力ファイルが存在するか確認
                if not os.path.exists(output_file):
                    raise RuntimeError("Output file was not created")
                
                # 新しいタンパク質オブジェクトを作成
                # 実際には、新しい構造データを読み込む必要があります
                prepared_protein = Protein(
                    id=protein.id,
                    structure=protein.structure,  # 本来は変換後の構造を読み込む
                    format=MoleculeFormat(type=FormatType.PDBQT),
                    properties=protein.properties.copy(),
                    path=output_file,  # 一時ファイルのパス（後で永続化が必要）
                    metadata=protein.metadata.copy(),
                    is_prepared=True,
                    preparation_method=prep_method,
                    chains=protein.chains.copy(),
                    has_water=False,  # 水分子は除去される
                    has_hydrogens=True,  # 水素原子は追加される
                    active_site_residues=protein.active_site_residues.copy()
                )
                
                return prepared_protein
                
        except Exception as e:
            logger.error(f"Error preparing receptor: {e}")
            raise RuntimeError(f"Error preparing receptor: {str(e)}")
    
    def calculate_properties(self, compound: Compound) -> List[MoleculeProperty]:
        """化合物の物理化学的特性を計算する
        
        Args:
            compound: 特性を計算する化合物
            
        Returns:
            計算された特性のリスト
            
        Raises:
            ValueError: 計算に失敗した場合
        """
        # この実装では、簡易的な特性計算のみを行う
        # 実際には、OpenBabelやRDKitなどのライブラリを使って計算する
        
        # プレースホルダーとしての実装
        properties = []
        
        # 分子量（ダミー値）
        from docking_automation.domain.molecule.value_object.molecule_property import MoleculeProperty, PropertyType
        properties.append(MoleculeProperty(
            name="molecular_weight",
            value=500.0,
            property_type=PropertyType.MOLECULAR_WEIGHT,
            unit="g/mol"
        ))
        
        return properties
    
    
    def add_hydrogens(self, molecule: MoleculeUnion, ph: float = 7.4) -> MoleculeUnion:
        """分子に水素原子を追加する
        
        Args:
            molecule: 水素原子を追加する分子
            ph: 溶液のpH値（デフォルト: 7.4）
            
        Returns:
            水素原子が追加された分子
            
        Raises:
            ValueError: 無効なパラメータが指定された場合
            RuntimeError: 処理に失敗した場合
        """
        # この実装では、OpenBabelを使って水素を追加する
        
        try:
            # 一時ディレクトリの作成
            with tempfile.TemporaryDirectory() as temp_dir:
                # 入力ファイルの準備
                input_file = molecule.path
                if not input_file or not os.path.exists(input_file):
                    raise ValueError("Input file is required")
                
                # 出力ファイルパスの設定
                output_file = os.path.join(temp_dir, f"{molecule.id}_h.{molecule.format.type.get_extension()}")
                
                # OpenBabelコマンドを構築
                cmd = [
                    self.obabel_path,
                    input_file,
                    "-O", output_file,
                    "-h",  # 水素原子を追加
                    f"-p {ph}"  # pH値を指定
                ]
                
                # コマンドを実行
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False
                )
                
                # 結果を確認
                if result.returncode != 0:
                    logger.error(f"Adding hydrogens failed: {result.stderr}")
                    raise RuntimeError(f"Adding hydrogens failed: {result.stderr}")
                
                # 出力ファイルが存在するか確認
                if not os.path.exists(output_file):
                    raise RuntimeError("Output file was not created")
                
                # 新しい分子オブジェクトを作成
                # 型に基づいて適切なオブジェクトを作成
                if isinstance(molecule, Compound):
                    # Compound型を返す
                    compound_result = Compound(
                        id=molecule.id,
                        structure=molecule.structure,  # 本来は変換後の構造を読み込む
                        format=molecule.format,
                        properties=molecule.properties.copy(),
                        path=output_file,
                        metadata=molecule.metadata.copy(),
                        is_prepared=molecule.is_prepared,
                        preparation_method=molecule.preparation_method
                    )
                    return compound_result
                elif isinstance(molecule, Protein):
                    # Protein型を返す（明示的に型を指定）
                    protein_result: Protein = Protein(
                        id=molecule.id,
                        structure=molecule.structure,  # 本来は変換後の構造を読み込む
                        format=molecule.format,
                        properties=molecule.properties.copy(),
                        path=output_file,
                        metadata=molecule.metadata.copy(),
                        is_prepared=molecule.is_prepared,
                        preparation_method=molecule.preparation_method,
                        chains=molecule.chains.copy(),
                        has_water=molecule.has_water,
                        has_hydrogens=True,  # 水素原子が追加された
                        active_site_residues=molecule.active_site_residues.copy()
                    )
                    return protein_result
                else:
                    raise TypeError(f"Unsupported molecule type: {type(molecule)}")
                
        except Exception as e:
            logger.error(f"Error adding hydrogens: {e}")
            raise RuntimeError(f"Error adding hydrogens: {str(e)}")
    
    def remove_water(self, protein: Protein) -> Protein:
        """タンパク質から水分子を除去する
        
        Args:
            protein: 水分子を除去するタンパク質
            
        Returns:
            水分子が除去されたタンパク質
            
        Raises:
            RuntimeError: 処理に失敗した場合
        """
        try:
            # 一時ディレクトリの作成
            with tempfile.TemporaryDirectory() as temp_dir:
                # 入力ファイルの準備
                input_file = protein.path
                if not input_file or not os.path.exists(input_file):
                    raise ValueError("Input file is required")
                
                # 出力ファイルパスの設定
                output_file = os.path.join(temp_dir, f"{protein.id}_nowater.{protein.format.type.get_extension()}")
                
                # 水分子を除去
                if not self._remove_water(input_file, output_file):
                    raise RuntimeError("Failed to remove water molecules")
                
                # 新しいタンパク質オブジェクトを作成
                result_protein = Protein(
                    id=protein.id,
                    structure=protein.structure,  # 本来は変換後の構造を読み込む
                    format=protein.format,
                    properties=protein.properties.copy(),
                    path=output_file,
                    metadata=protein.metadata.copy(),
                    is_prepared=protein.is_prepared,
                    preparation_method=protein.preparation_method,
                    chains=protein.chains.copy(),
                    has_water=False,  # 水分子が除去された
                    has_hydrogens=protein.has_hydrogens,
                    active_site_residues=protein.active_site_residues.copy()
                )
                
                return result_protein
                
        except Exception as e:
            logger.error(f"Error removing water: {e}")
            raise RuntimeError(f"Error removing water: {str(e)}")
    
    def get_supported_methods(self) -> Dict[str, str]:
        """このサービスがサポートする準備方法のリストを取得する
        
        Returns:
            サポートされている準備方法の辞書 {method_name: description}
        """
        return {
            "gasteiger": "Gasteiger電荷を使用した準備",
            "default": "デフォルトの準備手順",
        }
    
    def _convert_to_pdb(self, input_file: str, output_file: str) -> bool:
        """ファイルをPDB形式に変換する
        
        Args:
            input_file: 入力ファイルのパス
            output_file: 出力ファイルのパス
            
        Returns:
            変換が成功したかどうか
        """
        try:
            # OpenBabelコマンドを構築
            cmd = [
                self.obabel_path,
                input_file,
                "-O", output_file,
                "-h"  # 水素原子を追加
            ]
            
            # コマンドを実行
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )
            
            # 結果を確認
            if result.returncode != 0:
                logger.error(f"Conversion to PDB failed: {result.stderr}")
                return False
            
            # 出力ファイルが存在するか確認
            return os.path.exists(output_file)
            
        except Exception as e:
            logger.error(f"Error converting to PDB: {e}")
            return False
    
    def _remove_water(self, input_file: str, output_file: str) -> bool:
        """PDBファイルから水分子を除去する
        
        Args:
            input_file: 入力ファイルのパス
            output_file: 出力ファイルのパス
            
        Returns:
            処理が成功したかどうか
        """
        try:
            # 水分子を含まない行のみを抽出
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    # HETATM行のうち、HOHまたはWATを含まないものと、それ以外の行を抽出
                    if not (line.startswith("HETATM") and ("HOH" in line or "WAT" in line)):
                        outfile.write(line)
            
            return True
            
        except Exception as e:
            logger.error(f"Error removing water: {e}")
            return False