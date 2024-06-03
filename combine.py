import subprocess
import sys
import os

njob=int(sys.argv[1])   #352

cmd="> result/all_result.pdbqt"
subprocess.run(cmd, shell=True)

for i in range(njob):
    if os.path.exists("result/result"+str(i+1)+".pdbqt") != 1:
        continue
    cmd="cat result/result"+str(i+1)+".pdbqt >> result/all_result.pdbqt"
    subprocess.run(cmd, shell=True)

