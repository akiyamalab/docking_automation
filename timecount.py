#from datetime import datetime
import datetime
import sys

inputfile=open(sys.argv[1],'r')
lines=inputfile.readlines()
start_list=[]
end_list=[]


for line in lines:
    u_line=line.split()
    if u_line[1][0:5] == '53355':
        start_list.append(u_line[5].split('T')[0]+' '+u_line[5].split('T')[1])
        end_list.append(u_line[6].split('T')[0]+' '+u_line[6].split('T')[1])
start=sorted(start_list)[0]
end=sorted(end_list,reverse=True)[0]

starttime = datetime.datetime.strptime(start, '%Y-%m-%d %H:%M:%S')
endtime = datetime.datetime.strptime(end, '%Y-%m-%d %H:%M:%S')

ed_st=endtime-starttime
ed_st=ed_st.total_seconds()

alltime=0
for i in range(352):
    slurm_file=open('slurm-53355_'+str(i+1)+'.out','r')
    slines=slurm_file.readlines()
    for sline in slines:
        if sline[0:7] == 'alltime':
            alltime=alltime+int(sline.split()[2])

print(ed_st)
print("TIME= "+str((alltime/352)/ed_st))

