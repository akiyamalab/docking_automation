#manual

・ファイルは全て同じディレクトリ以下にあること
・vinaやprepare_ligandやprepare_receptorやslurmはインストールされていること（パスはauto_docking.py内で書き換えてください）

・実際にコマンドを打って実行するプログラムはdevide_sdf.pyとslurm.sh

#devide_sdf.py（sdfファイルと分割数を入力として分割数に応じたsdfファイルを作成するプログラム）
　入力：第一引数　sdfファイル
　　　　第二引数　分割数（500並列で実行するときは500を指定）
　出力：分割数がNの場合、カレントディレクトリ内に、PubChem_c_1.sdf〜PubChem_c_N.sdfのファイルが作成される


#slurm.sh（キューイングシステムにより、分割数に応じてドッキングを並列処理で実行するプログラム）
　入力：引数はないが、分割数と入力sdfファイルに応じてソースコード内の引数を書き換える必要がある
　　　　・分割数がNの時、３行目は　#SBATCH --array 1-N
　　　　・７行目のcount.pyの引数は　第一引数が分割前のsdfファイル、第二引数が分割数N
　出力：分割数Nの時、プログラム内で自動で作成されるresultディレクトリ内に、
　　　　a-result1.pdbqt〜a-resultN.pdbqtが作成される
