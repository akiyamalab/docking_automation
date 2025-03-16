from docking_automation.molecule import Protein
from docking_automation.molecule import CompoundSet
from docking_automation.docking import AutoDockVina, GridBox
from docking_automation.docking.autodockvina_docking import AutoDockVinaParameters

def main():
    """
    シンプルなドッキング計算の例
    """
    # 1. 分子データの読み込み（IDは自動生成）
    protein = Protein("path/to/protein.pdb")
    compound_set = CompoundSet("path/to/compounds.sdf")
    
    # 2. グリッドボックスの定義（タプルを使用）
    grid_box = GridBox(center=(10.0, 10.0, 10.0), size=(20.0, 20.0, 20.0))
    
    # 3. ドッキング計算（前処理は内部で自動実行）
    docking_tool = AutoDockVina()
    # AutoDockVinaParametersを作成して渡す
    additional_params = AutoDockVinaParameters()
    results = docking_tool.run_docking(protein, compound_set, grid_box, additional_params)
    
    # 4. 結果の取得と解析
    top_hits = results.get_top(10)
    
    # 結果の表示
    print(f"Top {len(top_hits)} hits:")
    for i, result in enumerate(top_hits):
        print(f"{i+1}. Score: {result.get_score()}")

if __name__ == "__main__":
    main()