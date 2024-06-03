import subprocess
import sys
import os

njob=int(sys.argv[1])   #352

cmd="> all_result_taioudata.txt"
subprocess.run(cmd, shell=True)

for i in range(njob):
    if os.path.exists('njob_taiou_'+str(i+1)+'.txt') != 1:
        continue
    cmd="cat njob_taiou_"+str(i+1)+".txt >> all_result_taioudata.txt"
    subprocess.run(cmd, shell=True)
