手順１、２は352並列で実行（10万件であれば10時間程度で終了する）

1,   作業ディレククトリ下のpre_shell_sbatch.shを実行することでsdfファイルからpdbqtファイルを複数生成
　　（コード中のinput_100000.sdfの部分を入力のsdfファイル名に書き換える）

2,   1のジョブが全て終了した後に、作業ディレクトリ下のshell_sbatch.shを実行することでドッキング計算及び計算時間の計測
　　を行う

3,   2のジョブが全て終了した後に作業ディレクトリ上で
　　　　python combine.py 352
             python combine_taiou.py 352
　　を実行すると、
　　作業ディレクトリ下のresultディレクトリ下にall_result.pdbqt（ドッキング結果を格納）と、
       作業ディレクトリ下にall_result_taioudata.txt（入力情報と計算時間を記録）が生成される。
　　
　　<all_result_taioudata.txt>
         [化合物識別タグ],[原子数],[回転可能な結合数],[分子量],[ドッキング計算時間]で構成される行が入力化合物の数の分ある

4,   作業ディレクトリ下のresultディレクトリに移動し、3で生成したall_result.pdbqtを引数として
　　　　python    sort_pdbqt.py    all_result.pdbqt
　　を実行すると、resultディレクトリ下にsorted-all_result.pdbqtが生成される。（スコアが良い順にソート）

5,   4で生成したsorted-all_resultを引数として、resultディレクトリ下で、
　　　　obabel  sorted-all_result.pdbqt  -O  sorted-all_result.sdf
　　を実行し、sdfファイルに変換（resultディレクトリ下に生成）

6,   5で生成したsorted-all_result.sdfに対しresultディレクトリ下で、
　　　　python  top_N_pull.py   sorted-all_result.sdf   N
      を実行し、resultディレクトリ下に、上位N件だけ抜き出したsdfファイルを生成

7,   作業ディレクトリに戻り、
             sacct -S 2019-04-01 -o User,JobID,Partition,NNodes,Submit,Start,End,Elapsed,State -X > alllog.txt
　　を実行した後、
　　　　python   timecount.py   alllog.txt
　　を実行すると
　　　　1行目に、初めのジョブが開始してから最後のジョブが終わるまでの時間（秒）（並列合計時間）と、
　　　　2行目に、（直列合計時間/352）/並列合計時間 　　が表示される（効率）
