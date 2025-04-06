#!/usr/bin/env python3
"""
タンパク質セグメンテーションスクリプト

単一のタンパク質を入力として、ProteinSegmentationServiceを使って
セグメント化された複数のタンパク質構造を出力します。
"""

import argparse
import os
from pathlib import Path
from typing import Any, Dict

from docking_automation.domain.services.protein_segmentation_service import (
    ProteinSegmentationService,
)
from docking_automation.molecule.protein import Protein


def main() -> int:
    """
    メイン関数

    Returns:
        int: 終了コード（0: 成功, 1: 失敗）
    """
    # コマンドライン引数の解析
    parser = argparse.ArgumentParser(description="タンパク質セグメンテーションツール")
    parser.add_argument("input_pdb", help="入力PDBファイルのパス")
    parser.add_argument("output_dir", help="出力ディレクトリのパス")
    parser.add_argument("--loop_min", type=int, default=20, help="最小ループ長")
    parser.add_argument("--helix_min", type=int, default=30, help="最小ヘリックス長")
    parser.add_argument("--fragment_min", type=int, default=5, help="最小フラグメント長")
    parser.add_argument("--domain_min", type=int, default=50, help="最小ドメイン長")
    parser.add_argument("--pLDDT_min", type=float, default=0, help="最小pLDDT値")
    parser.add_argument("--local_contact_range", type=int, default=5, help="局所接触範囲")
    args = parser.parse_args()

    # パスの設定
    input_path = Path(args.input_pdb)
    output_dir = Path(args.output_dir)

    # 入力ファイルの存在確認
    if not input_path.exists():
        print(f"エラー: 入力ファイル '{input_path}' が見つかりません。")
        return 1

    # 出力ディレクトリの作成
    os.makedirs(output_dir, exist_ok=True)

    # Proteinオブジェクトの作成
    protein = Protein.create(path=input_path)
    print(f"入力タンパク質: ID={protein.id}, パス={protein.path}")

    # AlphaCutterのオプション設定
    options: Dict[str, Any] = {
        "loop_min": args.loop_min,
        "helix_min": args.helix_min,
        "fragment_min": args.fragment_min,
        "domain_min": args.domain_min,
        "pLDDT_min": args.pLDDT_min,
        "local_contact_range": args.local_contact_range,
        "domain_out": True,
        "single_out": True,
    }

    # ProteinSegmentationServiceの作成と実行
    service = ProteinSegmentationService()
    print(f"タンパク質セグメンテーションを実行中...")
    segmented_proteins = service.segment(protein, options, output_dir)

    # 結果の表示
    print(f"\nセグメンテーション完了: {len(segmented_proteins)}個のセグメントが生成されました。")
    for i, seg_protein in enumerate(segmented_proteins, 1):
        print(f"セグメント {i}: ID={seg_protein.id}, パス={seg_protein.path}")

    # サマリーファイルの確認
    summary_file = output_dir / "AFCT-OUT_summary.csv"
    if summary_file.exists():
        print(f"\nサマリーファイル: {summary_file}")

    print(f"\n全ての出力ファイルは '{output_dir}' に保存されました。")
    return 0


if __name__ == "__main__":
    exit(main())
