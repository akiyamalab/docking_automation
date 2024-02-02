import sys
input_file_0=open(sys.argv[1],'r')
lines=input_file_0.readlines()
p=0
for line in lines:
    if line[0] == '$':
        p=p+1
print(p)
