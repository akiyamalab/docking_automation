#sys.argv[1]=sdf
#sys.argv[2]=n
#sys.argv[3]=h
import sys
import subprocess
h=int(sys.argv[2])
input_list=[]

import re
pattern=r'<(.*?)>'
pattern2=re.compile(r'[^a-zA-Z0-9_]')


input_file_0=open(sys.argv[1],'r',encoding="unicode_escape")#,encoding='utf-32')#,encoding="utf-8",errors='ignore')
lines=input_file_0.readlines()
p=0
for line in lines:
    if line[0] == '$':
        p=p+1

n=int(p/h)

input_file=open(sys.argv[1],'r',encoding="unicode_escape")#,encoding='utf-32')#,encoding="utf-8",errors='ignore')

for i in range(h):
    input_list.append(open('PubChem_c_'+str(i+1)+'.sdf','w'))#,encoding="unicode_escape"))#,encoding='utf-32'))
    input_list[i].truncate(0)


dc=0
fal=0

mk_flag=0

for line in input_file:
    if dc == 0:
        print(line)
    if line[0]=='\n' and fal==1:
        break
        
    fal=0
    if '<' in line:
        match=re.search(pattern,line)
        #print(match.group(1))
        #print(len(match.group(1)))
        if len(match.group(1)) == 4 and match.group(1) != 'CdId':
           print('aaaaaaaaa') 
           mk_flag=1
        else:
            input_list[int(dc/(n+1))].write(line)
    elif mk_flag == 1:
        mk_flag=0
    else:
        input_list[int(dc/(n+1))].write(line)
    if line[0] == '$':
        if dc != 0:
            pass
            #input_list[int(dc/(n+1))].write('\n')
        dc=dc+1
        fal=1
input_file.close()
input_file_0.close()
for i in range(h):
    input_list[i].close()

for i in range(h):
    cmd='python ccou.py pre_PubChem_c_'+str(i+1)+'.sdf '+str(i+1)
    #subprocess.run(cmd, shell=True)
    cmd='rm pre_PubChem_c_'+str(i+1)+'.sdf'
   # subprocess.run(cmd, shell=True)
