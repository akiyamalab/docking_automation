import subprocess
import sys
import os
input_file_0=open("result/"+sys.argv[1],'r')
if os.path.exists("result/"+sys.argv[1]) != 1:
    sys.exit(0)

taiou_file=open("taiou.txt",'r')  #全体の対応表
taioulist_file=open('taioulist'+sys.argv[2]+'.txt','r')

cmd = "touch "+"result/a-"+sys.argv[1]
subprocess.run(cmd, shell=True)
input_file_1=open("result/a-"+sys.argv[1],'w')
lines=input_file_0.readlines()
taioulines=taiou_file.readlines()


taioulist=taioulist_file.readlines()[0].split()

p=0
for line in lines:
    if line[0:5] == 'MODEL':
        kai=taioulist[p]
        for tline in taioulines:
            if tline.split(':')[0] == kai:
                 input_file_1.write("MODEL "+tline.split(':')[1].replace('\n','')+"\n")
        p=p+1
    else:
        input_file_1.write(line)

cmd = "rm result/"+sys.argv[1]
subprocess.run(cmd, shell=True)

input_file_0.close()
input_file_1.close()
taiou_file.close()
taioulist_file.close()

