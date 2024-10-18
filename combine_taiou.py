# 第二引数としてプロテイン数
import subprocess
import sys
import os

njob=int(sys.argv[1])   #352

npro=int(sys.argv[2])   #Number of protein

for u in range(npro):
    cmd="> all_result_taioudata_p"+str(u+1)+".txt"
    subprocess.run(cmd, shell=True)

for i in range(njob):
    for y in range(npro):
        if os.path.exists('njob_taiou_'+str(i+1)+'_p'+str(y+1)+'.txt') != 1:
            continue
        cmd="cat njob_taiou_"+str(i+1)+'_p'+str(y+1)+".txt >> all_result_taioudata_p"+str(y+1)+".txt"
        subprocess.run(cmd, shell=True)
        cmd="rm njob_taiou_"+str(i+1)+'_p'+str(y+1)+".txt"
        subprocess.run(cmd, shell=True)

