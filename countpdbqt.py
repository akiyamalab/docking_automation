import subprocess
import sys
import os
input_file_0=open(sys.argv[1],'r')
if os.path.exists(sys.argv[1]) != 1:
    sys.exit(0)



lines=input_file_0.readlines()




p=0
i=0
for line in lines:
    if line[0:8] == 'MODEL IT':
        p=p+1

    i=i+1


print(p)
print(i)

input_file_0.close()

