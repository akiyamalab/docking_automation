# ドッキング自動化パッケージ チュートリアル

このディレクトリには、docking_automationパッケージの使用方法を示す例示コードが含まれています。各サンプルスクリプトは、パッケージの異なる機能を紹介しています。

## 含まれるサンプルコード

1. **run_docking_example.py**
   - 新しいインポート構造を使用した基本的なドッキング計算の例
   - メモリ上でのモック分子とタンパク質の作成方法
   - ドッキングタスクのセットアップと実行方法

2. **run_simple_docking.py**
   - RDKitを使った単純な分子の生成とドッキング計算
   - 実際のSDFファイルからの分子読み込み
   - 結果の解析と要約の生成

3. **run_simple_docking_new.py**
   - 新しいインポート構造を使用した単純なドッキング計算の例
   - RDKitとMeekoを使用した分子生成と準備
   - docking_automationの構造体を使用した実装

4. **run_vina_docking.py**
   - 実際のVinaコマンドを使用したドッキング計算
   - タンパク質とリガンドの準備
   - グリッドボックスとパラメータの設定

5. **run_vina_direct.py**
   - 既存のPDBQTファイルを使用して直接Vinaでドッキング計算
   - 結果ファイルの解析と要約の生成
   - コマンドライン実行の例

6. **run_vina_docking_with_real_files.py**
   - 実際のファイルを使用したVinaドッキング計算の例
   - VinaDockingServiceクラスの使用方法
   - 複数レセプター対応の機能を示す実装

7. **run_meeko_vina_docking.py**
   - Meekoを使用したリガンド準備とVina計算
   - 複数化合物の一括処理
   - 結果の解析とランキング

8. **run_docking_workflow.py**
   - 新しいインポート構造を使用したドッキングワークフローの例
   - 複数化合物の処理とバッチ処理
   - 結果のランキングと要約

## 実行方法

各スクリプトはプロジェクトのルートディレクトリから以下のように実行できます：

```bash
python tutorial/run_docking_example.py
python tutorial/run_simple_docking.py
python tutorial/run_simple_docking_new.py
python tutorial/run_vina_docking.py
python tutorial/run_vina_direct.py
python tutorial/run_vina_docking_with_real_files.py
python tutorial/run_meeko_vina_docking.py
python tutorial/run_docking_workflow.py
```

## 必要な依存関係

これらのサンプルを実行するには以下のパッケージが必要です：

- RDKit
- Meeko
- AutoDock Vina (またはADVina)
- NumPy
- Matplotlib (結果の可視化用)

依存関係は以下のコマンドでインストールできます：

```bash
pip install rdkit meeko numpy matplotlib
```

AutoDock Vinaは別途インストールが必要です。

## サンプルの実行順序

初めての方は、以下の順序でサンプルを実行することをお勧めします：

1. `run_docking_example.py` - 基本的な使用方法の理解
2. `run_simple_docking.py` または `run_simple_docking_new.py` - 実際の分子でのドッキング
3. `run_vina_direct.py` - 既存ファイルを使用した直接計算
4. `run_vina_docking_with_real_files.py` - VinaDockingServiceの使用方法
5. `run_meeko_vina_docking.py` - 複数分子の処理
6. `run_docking_workflow.py` - 完全なワークフローの実行