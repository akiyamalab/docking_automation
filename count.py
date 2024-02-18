"""
python count.py [SDF] [number of subfiles]
"""
import sys
input_file_0=open(sys.argv[1],'r')
lines=input_file_0.readlines()
n_cmpds=0 # number of compounds in the SDF
n_jobs=sys.argv[2]
for line in lines:
    if line[0] == '$':
        n_cmpds=n_cmpds+1
n_cmpds_per_job=int(n_cmpds/int(n_jobs))
print(n_cmpds_per_job)
#print(p)
