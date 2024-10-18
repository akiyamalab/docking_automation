import sys

input_file=open(sys.argv[1],'r')
linelist=input_file.readlines()
input_file.close()

i=0
lastlist=[]
for line in linelist:
    if line[0:5] == 'MODEL':
        conlist=[]
        conlist.append(0)
        conlist.append(line)
    elif line[0:18] == 'REMARK VINA RESULT':
        if i != 0:
            lastlist.append(conlist)
        #conlist=[]

        li=line.split(':')[1]
        #conlist.append(int(li.split(' ')[1]))
        #conlist.append(line)
        #print(li.split(' '))

        conlist[0]=conlist[0]+float(li.lstrip().split(' ')[0])
        conlist[1]=conlist[1]+line

    else:
        conlist[1]=conlist[1]+line
    i=i+1

lastlist.append(conlist)
lastlist=sorted(lastlist,key=lambda x: x[0])

output_file=open(sys.argv[1].split('/')[0]+'/sorted-'+sys.argv[1].split('/')[1],'w')

for line in lastlist:
    output_file.write(line[1])

output_file.close()
