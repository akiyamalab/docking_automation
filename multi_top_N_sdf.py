# 第一引数にプロテイン数Nを指定
# 第二引数に上位M件のMを指定

import numpy as np
import sys
import subprocess
import os

N=int(sys.argv[1])
M=int(sys.argv[2])

for i in range(N):
    cmd="python sort_pdbqt.py result"+str(i+1)+"/all_result.pdbqt"
    subprocess.run(cmd, shell=True)
    cmd="obabel result"+str(i+1)+"/sorted-all_result.pdbqt -O result"+str(i+1)+"/sorted-all_result.sdf"
    subprocess.run(cmd, shell=True)
    cmd="python top_N_pull.py result"+str(i+1)+"/sorted-all_result.sdf "+str(M)
    subprocess.run(cmd, shell=True)

