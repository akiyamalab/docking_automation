import os
import tempfile
import logging
from typing import List, Tuple, TypeVar

from docking_automation.molecule.entity.molecule import Molecule
from docking_automation.molecule.entity.compound import Compound
from docking_automation.molecule.entity.protein import Protein
from docking_automation.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.molecule.service.format_converter_service import FormatConverterService

# 型変数の定義
MoleculeType = TypeVar('MoleculeType', bound=Molecule)

# ロガーの設定
logger = logging.getLogger(__name__)


class OpenBabelFormatConverter(FormatConverterService[MoleculeType]):
    """OpenBabelを使用した分子形式変換サービスの実装"""
    
    # フォーマット間の変換マッピング
    _FORMAT_MAP = {
        FormatType.SDF: "sdf",
        FormatType.PDB: "pdb",
        FormatType.PDBQT: "pdbqt",
        FormatType.MOL2: "mol2",
        FormatType.SMILES: "smi",
    }
    
    # サポートされている変換のリスト（入力形式、出力形式）
    _SUPPORTED_CONVERSIONS = [
        (FormatType.SDF, FormatType.PDB),
        (FormatType.SDF, FormatType.PDBQT),
        (FormatType.PDB, FormatType.PDBQT),
        (FormatType.PDBQT, FormatType.PDB),
        (FormatType.MOL2, FormatType.PDB),
        (FormatType.MOL2, FormatType.PDBQT),
        (FormatType.SMILES, FormatType.PDB),
        (FormatType.SMILES, FormatType.SDF),
    ]
    
    def __init__(self, obabel_path: str = "obabel"):
        """コンストラクタ
        
        Args:
            obabel_path: OpenBabelの実行ファイルパス
        """
        self.obabel_path = obabel_path
    
    def convert(self, molecule: MoleculeType, target_format: FormatType) -> MoleculeType:
        """分子を指定された形式に変換する
        
        Args:
            molecule: 変換する分子
            target_format: 変換先の形式
            
        Returns:
            変換された分子（元の分子と同じ型）
            
        Raises:
            ValueError: 変換元と変換先の形式が互換性がない場合
            IOError: ファイル操作に失敗した場合
        """
        # 安全にmolecule.formatとformat.typeにアクセス
        if not hasattr(molecule, 'format'):
            raise ValueError(f"Molecule {molecule.id if hasattr(molecule, 'id') else 'unknown'} has no format attribute")
        
        # formatがNoneまたはtype属性がない場合を処理
        molecule_format = molecule.format
        if molecule_format is None or not hasattr(molecule_format, 'type'):
            # ファイルパスがある場合は、そこから形式を推測
            if hasattr(molecule, 'path') and molecule.path:
                source_format = FormatType.from_path(molecule.path)
                if source_format == FormatType.UNKNOWN:
                    raise ValueError(f"Cannot determine format for molecule from path: {molecule.path}")
            else:
                raise ValueError(f"Molecule format is not specified and cannot be determined from path")
        else:
            source_format = molecule_format.type
        
        # 変換が必要ない場合は元の分子をそのまま返す
        if source_format == target_format:
            return molecule
        
        # 変換がサポートされているか確認
        if not self.can_convert(source_format, target_format):
            raise ValueError(
                f"Conversion from {source_format} to {target_format} is not supported"
            )
        
        # 入力ファイルがない場合は一時ファイルを作成
        input_path = molecule.path if hasattr(molecule, 'path') else None
        input_created = False
        
        if not input_path or not os.path.exists(input_path):
            # 一時ファイルを作成
            with tempfile.NamedTemporaryFile(
                suffix=f".{self._FORMAT_MAP[source_format]}", delete=False
            ) as temp:
                input_path = temp.name
                input_created = True
                # 将来的には分子構造をファイルに書き出す実装が必要
                # ここではプレースホルダー
                temp.write(b"Placeholder for molecule data")
        
        try:
            # 出力用の一時ファイルを作成
            with tempfile.NamedTemporaryFile(
                suffix=f".{self._FORMAT_MAP[target_format]}", delete=False
            ) as output_file:
                output_path = output_file.name
            
            # 変換を実行
            success = self.convert_file(input_path, output_path, target_format)
            
            if not success:
                raise IOError(f"Failed to convert from {source_format} to {target_format}")
            
            # 変換結果を読み込む
            # 実際には、ファイルから分子データを読み込む実装が必要
            # ここでは、仮の実装として新しいオブジェクトを作成
            
            # 安全にMoleculeFormatインスタンスを作成
            new_format = MoleculeFormat(type=target_format)
            logger.debug(f"Created new format: {new_format} for target format: {target_format}")
            
            # 元の分子と同じ型のコンストラクタでインスタンスを作成
            # これにより、型チェッカーのMoleculeType型制約を満たす
            constructor = type(molecule)
            logger.debug(f"Using constructor: {constructor.__name__} for molecule type: {type(molecule).__name__}")
            
            # 共通の引数を準備（安全にアクセス）
            common_args = {
                "id": getattr(molecule, 'id', None) or str(id(molecule)),
                "path": output_path,
            }
            
            # 構造情報をコピー（存在する場合）
            if hasattr(molecule, 'structure'):
                common_args["structure"] = molecule.structure
            
            # 形式情報を設定
            common_args["format"] = new_format
            
            # プロパティをコピー（存在する場合）
            if hasattr(molecule, 'properties') and hasattr(molecule.properties, 'copy'):
                common_args["properties"] = molecule.properties.copy()
            elif hasattr(molecule, 'properties'):
                common_args["properties"] = molecule.properties
            
            # メタデータをコピー（存在する場合）
            if hasattr(molecule, 'metadata') and hasattr(molecule.metadata, 'copy'):
                common_args["metadata"] = molecule.metadata.copy()
            elif hasattr(molecule, 'metadata'):
                common_args["metadata"] = molecule.metadata
            
            # 型に応じて追加の引数を準備
            additional_args = {}
            
            if isinstance(molecule, Compound):
                # Compoundクラスの属性を安全に取得
                # 型チェッカーを抑制しつつ、属性を安全に取得
                is_prepared = getattr(molecule, 'is_prepared', False)  # type: ignore
                preparation_method = getattr(molecule, 'preparation_method', None)  # type: ignore
                
                additional_args["is_prepared"] = is_prepared
                if preparation_method is not None:
                    additional_args["preparation_method"] = preparation_method
                
            elif isinstance(molecule, Protein):
                # Proteinクラスの属性を安全に取得
                # 型チェッカーを抑制しつつ、属性を安全に取得
                is_prepared = getattr(molecule, 'is_prepared', False)  # type: ignore
                preparation_method = getattr(molecule, 'preparation_method', None)  # type: ignore
                chains = getattr(molecule, 'chains', set())  # type: ignore
                has_water = getattr(molecule, 'has_water', False)  # type: ignore
                has_hydrogens = getattr(molecule, 'has_hydrogens', False)  # type: ignore
                active_site_residues = getattr(molecule, 'active_site_residues', [])  # type: ignore
                
                additional_args["is_prepared"] = is_prepared
                if preparation_method is not None:
                    additional_args["preparation_method"] = preparation_method
                additional_args["chains"] = chains.copy() if hasattr(chains, 'copy') else chains
                additional_args["has_water"] = has_water
                additional_args["has_hydrogens"] = has_hydrogens
                additional_args["active_site_residues"] = active_site_residues.copy() if hasattr(active_site_residues, 'copy') else active_site_residues
                
            else:
                raise TypeError(f"Unsupported molecule type: {type(molecule)}")
                
            common_args.update(additional_args)
            
            try:
                # 型安全なインスタンス生成
                # TypeVarの境界がMoleculeなので、constructorの戻り値の型はMoleculeTypeになる
                
                # コンストラクタの引数を調査し、サポートしているものだけを使用
                import inspect
                constructor_sig = inspect.signature(constructor)
                valid_args = {}
                
                for param_name in constructor_sig.parameters:
                    if param_name in common_args:
                        valid_args[param_name] = common_args[param_name]
                
                logger.debug(f"Creating new molecule instance with valid args: {list(valid_args.keys())}")
                converted = constructor(**valid_args)  # type: ignore
                logger.info(f"Successfully converted molecule {getattr(molecule, 'id', 'unknown')} from {source_format} to {target_format}")
                return converted
            except Exception as e:
                logger.error(f"Error creating new molecule instance: {e}")
                raise ValueError(f"Failed to create new molecule instance: {e}")
            
        finally:
            # 一時ファイルを削除
            if input_created and os.path.exists(input_path):
                try:
                    os.unlink(input_path)
                    logger.debug(f"Removed temporary input file: {input_path}")
                except Exception as e:
                    logger.warning(f"Failed to remove temporary input file {input_path}: {e}")
    
    def convert_file(self, input_path: str, output_path: str, target_format: FormatType) -> bool:
        """ファイルを指定された形式に変換する
        
        Args:
            input_path: 入力ファイルのパス
            output_path: 出力ファイルのパス
            target_format: 変換先の形式
            
        Returns:
            変換が成功したかどうか
            
        Raises:
            ValueError: 入力ファイルが無効な場合
            IOError: ファイル操作に失敗した場合
        """
        import subprocess
        
        if not os.path.exists(input_path):
            raise ValueError(f"Input file not found: {input_path}")
        
        try:
            # 入力ファイルの形式を拡張子から推測
            input_format = FormatType.from_path(input_path)
            
            # UNKNOWNの場合はファイル名から推測できないため、エラーを出す
            if input_format == FormatType.UNKNOWN:
                raise ValueError(f"Cannot determine format for file: {input_path}")
            
            # 変換がサポートされているか確認
            if not self.can_convert(input_format, target_format):
                raise ValueError(
                    f"Conversion from {input_format} to {target_format} is not supported"
                )
        except Exception as e:
            logger.error(f"Error determining file format: {e}")
            raise ValueError(f"Failed to process input file: {input_path}: {str(e)}")
        
        # OpenBabelコマンドを構築
        cmd = [
            self.obabel_path,
            input_path,
            "-O",
            output_path,
            "-h"  # 水素原子を追加
        ]
        
        # 特別なオプション（例：PDBQT形式への変換時）
        if target_format == FormatType.PDBQT:
            cmd.extend(["-xr"])  # 回転可能な結合を処理
        
        try:
            # コマンドを実行
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False  # エラーがあっても例外をスローしない
            )
            
            # 結果を確認
            if result.returncode != 0:
                logger.error(f"OpenBabel conversion failed: {result.stderr}")
                return False
            
            # 出力ファイルが存在するか確認
            if not os.path.exists(output_path):
                logger.error("Output file was not created")
                return False
            
            return True
            
        except Exception as e:
            logger.error(f"Error during conversion: {e}")
            return False
    
    def get_supported_conversions(self) -> List[Tuple[FormatType, FormatType]]:
        """このサービスがサポートする変換のリストを取得する
        
        Returns:
            サポートされている変換のリスト（入力形式、出力形式）のタプル
        """
        return self._SUPPORTED_CONVERSIONS.copy()
    
    def can_convert(self, source_format: FormatType, target_format: FormatType) -> bool:
        """指定された変換がサポートされているかどうかを判定する
        
        Args:
            source_format: 変換元の形式
            target_format: 変換先の形式
            
        Returns:
            変換がサポートされているかどうか
        """
        # 同じ形式の場合は常にサポート
        if source_format == target_format:
            return True
        
        # サポートされている変換のリストで確認
        return (source_format, target_format) in self._SUPPORTED_CONVERSIONS