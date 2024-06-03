import sys
input_file=open(sys.argv[1],'r')
N=int(sys.argv[2]) #COUNT
#njob=int(sys.argv[3]) #352
lines=input_file.readlines()
output_file=open('top'+str(N)+'_'+sys.argv[1],'w')
p=0
for line in lines:
    output_file.write(line)
    if line[0] == '$':
        p=p+1
        if p == N:
            break
input_file.close()
output_file.close()
