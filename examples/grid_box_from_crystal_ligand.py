"""
結晶リガンド（crystal_ligand.mol2）からGridBoxを生成するサンプルコード
"""

import tempfile
from pathlib import Path

# OpenBabelのインポート
import openbabel.pybel as pybel

from docking_automation.docking.grid_box import GridBox
from docking_automation.molecule.compound_set import CompoundSet


def convert_mol2_to_sdf(mol2_path: Path) -> Path:
    """
    mol2ファイルをsdfファイルに変換する

    Args:
        mol2_path: mol2ファイルのパス

    Returns:
        変換後のsdfファイルのパス
    """
    # 一時ファイルを作成
    with tempfile.NamedTemporaryFile(delete=False, suffix=".sdf") as temp_file:
        temp_path = Path(temp_file.name)

    # mol2ファイルを読み込む
    mol = next(pybel.readfile("mol2", str(mol2_path)))

    # sdfファイルとして保存（overwrite=Trueを指定）
    mol.write("sdf", str(temp_path), overwrite=True)

    return temp_path


def main():
    """
    メイン関数
    """
    # 結晶リガンドのパス
    crystal_ligand_path = Path("examples/input/ALDR/crystal_ligand.mol2")
    # mol2ファイルをsdfファイルに変換
    sdf_path = convert_mol2_to_sdf(crystal_ligand_path)
    # CompoundSetの作成
    compound_set = CompoundSet(path=sdf_path)
    # GridBoxの生成
    grid_box = GridBox.from_compound(compound_set)

    # 結果の表示
    print("\n=== GridBox情報 ===")
    print(f"中心座標: {grid_box.center}")
    print(f"サイズ: {grid_box.size}")


if __name__ == "__main__":
    main()
