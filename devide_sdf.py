#sys.argv[1]=sdf
#sys.argv[2]=n
#sys.argv[3]=h
import sys
n_subfiles=int(sys.argv[2])
input_list=[]



input_file_0=open(sys.argv[1],'r')
lines=input_file_0.readlines()
n_cmpds=0
for line in lines:
    if line[0] == '$':
        n_cmpds=n_cmpds+1

n_cmpds_per_file=int(n_cmpds/n_subfiles)

input_file=open(sys.argv[1],'r')

for i in range(n_subfiles):
    input_list.append(open('PubChem_c_'+str(i+1)+'.sdf','a'))
    input_list[i].truncate(0)


#output_file1=open('PubChem_c_1.sdf','w')
#output_file2=open('PubChem_c_2.sdf','w')
cmpd_idx=0 # n-th compound, 0-origin
fal=0


for line in input_file:
    if line[0]=='\n' and fal==1:
        break
    fal=0
#    if dc <= 49:
#        output_file1.write(line)
#    else:
#        output_file2.write(line)

#    if int(dc/n)==h:
#            input_list[h-1].write(line)
#    else:
    input_list[int(cmpd_idx/(n_cmpds_per_file+1))].write(line)
    if line[0] == '$':
        cmpd_idx=cmpd_idx+1
        fal=1
input_file.close()
input_file_0.close()
for i in range(n_subfiles):
    input_list[i].close()
