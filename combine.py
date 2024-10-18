import subprocess
import sys
import os

njob=int(sys.argv[1])   #352

npro=int(sys.argv[2])   #Number of protein

for p in range(npro):
    cmd="> result"+str(p+1)+"/all_result.pdbqt"
subprocess.run(cmd, shell=True)

for i in range(njob):
    for p in range(npro):
        if os.path.exists("result"+str(p+1)+"/result"+str(i+1)+".pdbqt") != 1:
            continue
        cmd="cat result"+str(p+1)+"/result"+str(i+1)+".pdbqt >> result"+str(p+1)+"/all_result.pdbqt"
        subprocess.run(cmd, shell=True)

