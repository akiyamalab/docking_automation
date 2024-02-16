import sys
input_file_0=open(sys.argv[1],'r')
input_file_1=open('input_'+sys.argv[2]+'.sdf','w')
lines=input_file_0.readlines()
p=1
for line in lines:
    if p > int(sys.argv[2]):
        break
    input_file_1.write(line)
    if line[0] == '$':
        p=p+1
input_file_0.close()
input_file_1.close()
