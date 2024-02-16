#manual

・ファイルは全て同じディレクトリ以下にあること
・vinaやprepare_ligandやprepare_receptorやslurmはインストールされていること（パスはauto_docking.py内で書き換えてください）

・実際にコマンドを打って実行するプログラムはmake_v_num.pyとshell_sbatch.sh
＜手順＞
①make_v_num.pyを引数を指定して実行
②shell_sbatch.shを実行




#make_v_num.py（第一引数で指定したsdfファイルから第二引数で指定した数M件分だけ前から抜き出したsdfファイルを作成するプログラム）
　入力：第一引数　sdfファイル
　　　　第二引数　使用するデータ件数（10000化合物で実験するときは10000を指定）
　出力：指定したデータ件数がMの場合、カレントディレクトリ内に、input_M.sdfというファイルが作成される

#shell_sbatch.sh(内部でsdfからpdbへの変換や、sbatchによるジョブ投入を行うプログラム）
　入力：引数の指定は必要なし
　　　　入力に応じてコード内の引数を書き換える必要あり
　　　　・3行目のcount.pyの第一引数はmake_v_num.pyで作成したinput_M.sdfを指定
　　　　・6行目のconvert_sdf_to_pdb.pyの第一引数はmake_v_num.pyで作成したinput_M.sdfを指定
　　　　・3行目のcount.pyの第二引数は分割数（並列数）Nを指定
　　　　・6行目のconvert_sdf_to_pdb.pyの第三引数は分割数（並列数）Nを指定

#slurm.sh（キューイングシステムにより、分割数に応じてドッキングを並列処理で実行するプログラム）
*shell_sbatch.sh内で実行されるのでコマンドを打つ必要なし
　入力：分割数に応じてソースコード内の引数を書き換える必要がある
　　　　・分割数がNの時、３行目は　#SBATCH --array 1-N
　出力：分割数Nの時、プログラム内で自動で作成されるresultディレクトリ内に、
　　　　a-result1.pdbqt〜a-resultN.pdbqtが作成される
