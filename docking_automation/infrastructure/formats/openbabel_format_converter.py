import os
import tempfile
import logging
from typing import List, Tuple, TypeVar

from docking_automation.domain.molecule.entity.molecule import Molecule
from docking_automation.domain.molecule.entity.compound import Compound
from docking_automation.domain.molecule.entity.protein import Protein
from docking_automation.domain.molecule.value_object.molecule_format import MoleculeFormat, FormatType
from docking_automation.domain.molecule.service.format_converter_service import FormatConverterService

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
        # 変換が必要ない場合は元の分子をそのまま返す
        if molecule.format.type == target_format:
            return molecule
        
        # 変換がサポートされているか確認
        if not self.can_convert(molecule.format.type, target_format):
            raise ValueError(
                f"Conversion from {molecule.format.type} to {target_format} is not supported"
            )
        
        # 入力ファイルがない場合は一時ファイルを作成
        input_path = molecule.path
        input_created = False
        
        if not input_path or not os.path.exists(input_path):
            # 一時ファイルを作成
            with tempfile.NamedTemporaryFile(
                suffix=f".{self._FORMAT_MAP[molecule.format.type]}", delete=False
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
                raise IOError(f"Failed to convert from {molecule.format.type} to {target_format}")
            
            # 変換結果を読み込む
            # 実際には、ファイルから分子データを読み込む実装が必要
            # ここでは、仮の実装として新しいオブジェクトを作成
            
            new_format = MoleculeFormat(type=target_format)
            
            # 元の分子と同じ型のコンストラクタでインスタンスを作成
            # これにより、型チェッカーのMoleculeType型制約を満たす
            constructor = type(molecule)
            
            # 共通の引数を準備
            common_args = {
                "id": molecule.id,
                "structure": molecule.structure,  # 実際には変換後の構造を読み込む
                "format": new_format,
                "properties": molecule.properties.copy(),
                "path": output_path,
                "metadata": molecule.metadata.copy(),
            }
            
            # 型に応じて追加の引数を準備
            if isinstance(molecule, Compound):
                common_args.update({
                    "is_prepared": molecule.is_prepared,
                    "preparation_method": molecule.preparation_method
                })
            elif isinstance(molecule, Protein):
                common_args.update({
                    "is_prepared": molecule.is_prepared,
                    "preparation_method": molecule.preparation_method,
                    "chains": molecule.chains.copy(),
                    "has_water": molecule.has_water,
                    "has_hydrogens": molecule.has_hydrogens,
                    "active_site_residues": molecule.active_site_residues.copy()
                })
            else:
                raise TypeError(f"Unsupported molecule type: {type(molecule)}")
            
            # 型安全なインスタンス生成
            # TypeVarの境界がMoleculeなので、constructorの戻り値の型はMoleculeTypeになる
            converted = constructor(**common_args)  # type: ignore
            return converted
            
        finally:
            # 一時ファイルを削除
            if input_created and os.path.exists(input_path):
                os.unlink(input_path)
    
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
        
        # 入力ファイルの形式を拡張子から推測
        input_format = FormatType.from_path(input_path)
        
        # 変換がサポートされているか確認
        if not self.can_convert(input_format, target_format):
            raise ValueError(
                f"Conversion from {input_format} to {target_format} is not supported"
            )
        
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