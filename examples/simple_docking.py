import os
import sys
import argparse
from io import StringIO
from pathlib import Path
from docking_automation.molecule import Protein
from docking_automation.molecule import CompoundSet
from docking_automation.docking import AutoDockVina, GridBox
from docking_automation.docking.autodockvina_docking import AutoDockVinaParameters

# OpenBabelのwarningを完全に抑制するための設定
os.environ["BABEL_QUIET"] = "1"

def parse_args():
    """
    コマンドライン引数を解析する。
    
    Returns:
        argparse.Namespace: 解析されたコマンドライン引数
    """
    parser = argparse.ArgumentParser(description='ALDRタンパク質に対するドッキング計算を実行する')
    parser.add_argument('--decoys', action='store_true', help='デコイ化合物を使用する（デフォルトは活性化合物）')
    parser.add_argument('--exhaustiveness', type=int, default=4, help='探索の徹底度（デフォルト: 4）')
    parser.add_argument('--num-modes', type=int, default=3, help='出力するポーズの数（デフォルト: 3）')
    parser.add_argument('--energy-range', type=float, default=3.0, help='出力するポーズのエネルギー範囲（デフォルト: 3.0）')
    parser.add_argument('--top', type=int, default=10, help='表示する上位ヒット数（デフォルト: 10）')
    return parser.parse_args()

def main():
    """
    シンプルなドッキング計算の例
    ALDR（アルドース還元酵素）のデータセットを使用
    """
    # コマンドライン引数を解析
    args = parse_args()
    
    # 標準エラー出力をキャプチャするための設定
    # 元の標準エラー出力を保存
    original_stderr = sys.stderr
    # 標準エラー出力をStringIOにリダイレクト
    error_output = StringIO()
    sys.stderr = error_output
    
    try:
        # 現在のスクリプトのディレクトリを取得
        script_dir = Path(os.path.dirname(os.path.abspath(__file__)))
        
        # 1. 分子データの読み込み（IDは自動生成）
        # receptor.pdbファイルを使用する（正しいタンパク質ファイル）
        protein_path = script_dir / "input" / "ALDR" / "receptor.pdb"
        print(f"タンパク質ファイルの絶対パス: {protein_path.absolute()}")
        
        # 使用する化合物セットを選択（活性化合物またはデコイ化合物）
        use_actives = not args.decoys
        
        if use_actives:
            # 活性化合物を使用
            compound_path = script_dir / "input" / "ALDR" / "actives_final.sdf.gz"
            print("活性化合物を使用します。")
        else:
            # デコイ化合物を使用
            compound_path = script_dir / "input" / "ALDR" / "decoys_final.sdf.gz"
            print("デコイ化合物を使用します。")
        
        print(f"タンパク質ファイル: {protein_path}")
        print(f"化合物ファイル: {compound_path}")
        
        # Proteinオブジェクトの初期化（prepare_receptorが実行される部分）
        try:
            print("タンパク質の前処理を開始します...")
            protein = Protein(protein_path)
            print("タンパク質の前処理が完了しました")
        except Exception as e:
            print(f"タンパク質の前処理中にエラーが発生しました: {e}")
            # エラー出力を即座に表示
            error_content = error_output.getvalue()
            if error_content:
                print("\n=== prepare_receptorのエラー出力 ===")
                print(error_content)
                print("=== エラー出力終了 ===\n")
            raise  # エラーを再度発生させて処理を中断
            
        compound_set = CompoundSet(compound_path)
        
        # 化合物数を表示
        print(f"化合物数: {compound_set.get_compound_count()}")
        
        # 2. グリッドボックスの定義
        # 結晶構造のリガンド位置を中心とする
        # 実際のプロジェクトでは、結晶構造のリガンド位置や既知の活性部位情報から設定する
        grid_box = GridBox(center=(15.0, 23.0, 36.0), size=(20.0, 20.0, 20.0))
        
        # 3. ドッキング計算（前処理は内部で自動実行）
        docking_tool = AutoDockVina()
        
        # AutoDockVinaParametersを作成して渡す
        # コマンドライン引数から設定を取得
        additional_params = AutoDockVinaParameters(
            exhaustiveness=args.exhaustiveness,  # 探索の徹底度
            num_modes=args.num_modes,           # 出力するポーズの数
            energy_range=args.energy_range      # 出力するポーズのエネルギー範囲
        )
        
        print(f"ドッキングパラメータ: exhaustiveness={args.exhaustiveness}, num_modes={args.num_modes}, energy_range={args.energy_range}")
        print("ドッキング計算を開始します...")
        results = docking_tool.run_docking(protein, compound_set, grid_box, additional_params)
        
        # 4. 結果の取得と解析
        top_hits = results.get_top(args.top)
        
        # 結果の表示
        print(f"\nTop {len(top_hits)} hits:")
        for i, result in enumerate(top_hits):
            print(f"{i+1}. Score: {result.docking_score}, Compound: {result.compound_set_id}_{result.compound_index}")
            print(f"   Pose file: {result.result_path}")
    except Exception as e:
        print(f"エラーが発生しました: {e}")
    finally:
        # キャプチャしたエラー出力を表示
        error_content = error_output.getvalue()
        if error_content:
            print("\n=== 標準エラー出力 ===")
            print(error_content)
            print("=== エラー出力終了 ===\n")
        
        # 標準エラー出力を元に戻す
        sys.stderr = original_stderr

if __name__ == "__main__":
    main()