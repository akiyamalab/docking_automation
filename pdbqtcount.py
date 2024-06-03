import sys
from openbabel import openbabel
import os
import datetime
import sys
import csv
#import pybel


def get_pdbqt_files(directory):
    files = [f for f in os.listdir(directory) if f.startswith('output') and f.endswith('.pdbqt')]
    return files

pdbqt_files = get_pdbqt_files('.')
print(len(pdbqt_files))

for filename in pdbqt_files:
    m=0
    input_file=open(filename,'r')
    lines=input_file.readlines()
    for line in lines:
        if line[0:12]=='REMARK  Name':
              if m == 1:
                 print(filename) 
              m=m+1
