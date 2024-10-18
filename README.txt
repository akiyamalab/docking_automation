＊前提とする入力によって手順0から実行するか１から実行するかは適宜変更
＜手順０から始める際の前提＞
テストデータであるoriginalフォルダのような形式（中に複数のプロテイン名のフォルダと、各フォルダ内にもpdbとmol2ファイルを持つ）のフォルダがあること
input_proteinディレクトリ、input_confディレクトリは用意する必要なし
＊confファイルにおけるタンパク質中のグリッドの中心座標は、タンパク質残基座標の中心にしてある

＜手順１から始める際の前提＞
input_proteinディレクトリにプロテインファイルを配置・ファイル名はprotein1.pdb〜proteinN.pdb
Input_confディレクトリにグリッドファイルを配置・ファイル名はプロテインファイルと対応するようにautodock1.conf〜autodockN.conf


手順１、２は352並列で実行（10万件であれば10時間程度で終了する）

＊手順４〜６は各タンパク質に対して個別に行うが、
　　　　multip_top_N_sdf.py  （プロテインの数）　（M）
　　を実行することで、すべてのタンパク質に対して４〜６の手順が実行される。（Mは上位M件をsdfとして抽出したい時）


0,   テストデータであるoriginalフォルダのような形式（中に複数のプロテイン名のフォルダと、各フォルダ内にもpdbとmol2ファイルを持つ）のフォルダのパスを指定し、
　　　　python make_pro_conf.py 　（パス）
　　を実行すると、作業ディレクトリ下に各プロテイン名と１〜Nの対応関係を記録したpro_number_log.txtが作成される。
　　また、上記の対応関係に従ってprotein1.pdb〜proteinN.pdbがinput_proteinディレクトリに、autodock1.conf〜autodockN.confがinput_confディレクトリに作成される。

1,   作業ディレククトリ下のpre_shell_sbatch.shを実行することでsdfファイルからpdbqtファイルを複数生成
　　（コード中のinput_100000.sdfの部分を入力のsdfファイル名に書き換える）

2,   1のジョブが全て終了した後に、作業ディレクトリ下のshell_sbatch.shを実行することでドッキング計算及び計算時間の計測
　　を行う

3,   2のジョブが全て終了した後に作業ディレクトリ上で
　　　　python combine.py 352  (proteinの数)
             python combine_taiou.py 352  (proteinの数)
　　を実行すると、
　　作業ディレクトリ下の各タンパク質のresultディレクトリ下にall_result.pdbqt（ドッキング結果を格納）と、
       作業ディレクトリ下にall_result_taioudata_p1.txt~all_result_taioudata_pN.txt（入力情報と計算時間を記録・Nはタンパク質の数）が生成される。
　　
　　<all_result_taioudata.txt>
         [化合物識別タグ],[原子数],[回転可能な結合数],[分子量],[ドッキング計算時間]で構成される行が入力化合物の数の分ある

4,   作業ディレクトリ上で、3で生成した対応するタンパク質Nの、resultN/all_result.pdbqtを引数として
　　　　python    sort_pdbqt.py    resultN/all_result.pdbqt
　　を実行すると、対応するタンパク質のresultディレクトリ下にsorted-all_result.pdbqtが生成される。（スコアが良い順にソート）

5,   4で生成したsorted-all_resultを引数として、対応するタンパク質のresultディレクトリ下で、
　　　　obabel  resultN/sorted-all_result.pdbqt  -O  resultN/sorted-all_result.sdf
　　を実行し、sdfファイルに変換（対応するタンパク質のresultディレクトリ下に生成）

6,   5で生成したsorted-all_result.sdfに対し、作業ディレクトリ下で、
　　　　python  top_N_pull.py   resultN/sorted-all_result.sdf   N
      を実行し、対応するタンパク質のresultディレクトリ下に、上位N件だけ抜き出したsdfファイルを生成

7,   作業ディレクトリに戻り、
             sacct -S 2019-04-01 -o User,JobID,Partition,NNodes,Submit,Start,End,Elapsed,State -X > alllog.txt
　　を実行した後、
　　　　python   timecount.py   alllog.txt
　　を実行すると
　　　　1行目に、初めのジョブが開始してから最後のジョブが終わるまでの時間（秒）（並列合計時間）と、
　　　　2行目に、（直列合計時間/352）/並列合計時間 　　が表示される（効率）


      ＊timecount.pyのコード中でslurm IDを指定している部分が２箇所あるのでそこは適宜変更の必要あり
