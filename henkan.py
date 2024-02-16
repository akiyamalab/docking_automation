import subprocess
import sys
import os
input_file_0=open("result/"+sys.argv[1],'r')
if os.path.exists("result/"+sys.argv[1]) != 1:
    sys.exit(0)

cmd = "touch "+"result/a-"+sys.argv[1]
subprocess.run(cmd, shell=True)
input_file_1=open("result/a-"+sys.argv[1],'w')
lines=input_file_0.readlines()
p=1
for line in lines:
    if line[0:5] == 'MODEL':
        input_file_1.write("MODEL "+str(p)+"\n")
    else:
        input_file_1.write(line)
    p=p+1
cmd = "rm result/"+sys.argv[1]
subprocess.run(cmd, shell=True)
