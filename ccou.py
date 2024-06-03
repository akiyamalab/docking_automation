import sys

input_file_0=open(sys.argv[1],'r',encoding="unicode_escape")
lines=input_file_0.readlines()
p=0
h=sys.argv[2]
outp=open('PubChem_c_'+h+'.sdf','w')
flag=0
for line in lines:
    if '<>' in line:
        flag=1
    elif flag == 1:
        flag=0
    else:
        outp.write(line)

