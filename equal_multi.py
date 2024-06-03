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



def get_molecule_weight(pdbqt_file):
    # Open BabelのMoleculeオブジェクトを作成
    mol = openbabel.OBMol()

    # pdbqtファイルを読み込む
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("pdbqt", "pdbqt")
    conv.ReadFile(mol, pdbqt_file)

    # 分子量を取得
    mol_wight = mol.GetMolWt()

    return mol_wight


def get_molecule_length(pdbqt_file):
    # Open BabelのMoleculeオブジェクトを作成
    mol = openbabel.OBMol()

    # pdbqtファイルを読み込む
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("pdbqt", "pdbqt")
    conv.ReadFile(mol, pdbqt_file)

    # 分子の原子の数を取得
    num_atoms = mol.NumAtoms()

    return num_atoms

def count_active_torsions(file_path):
    obConversion = openbabel.OBConversion()
    obConversion.SetInFormat("pdbqt")

    mol = openbabel.OBMol()
    obConversion.ReadFile(mol, file_path)

    torsion_count = 0
    for bond in openbabel.OBMolBondIter(mol):
        if bond.IsRotor():
            torsion_count += 1

    return torsion_count

current_directory = os.getcwd()  # 現在のディレクトリを取得
pdbqt_files = get_pdbqt_files(current_directory)

output_over_file=open('over.txt','w')
output_taiou_file=open('data_taiou.txt','w')

sum_len=0
num=0

for filename in pdbqt_files:
    input_file=open(filename,'r')
    lines=input_file.readlines()
    for line in lines:
        if line[0:12]=='REMARK  Name':
            num=num+1
            #print(line.split())
            title=line.split()[3]
            mol_len=int(get_molecule_length(filename))
            sum_len=sum_len+mol_len
            mol_tor=int(count_active_torsions(filename))
            output_taiou_file.write(title+','+str(mol_len)+','+str(mol_tor)+'\n')

            weight=get_molecule_weight(filename)
            if float(weight) > 700:
                output_over_file.write(title+'\n')

output_over_file.close()
output_taiou_file.close()

#avr_len=int(sum_len/num)+1
avr_len=int(sum_len/352)+1

#print('sum_len='+str(sum_len))
#print('avr_len='+str(avr_len))

output_job_file=[]
for j in range(352):
    output_job_file.append(open('joblist_'+str(j+1)+'.txt','w'))

i=0
comtemp_len=0
t_sum_len=0
lom=0
for filename in pdbqt_files:
    output_job_file[lom].write(filename+'\n')
    mol_len=int(get_molecule_length(filename))
    comtemp_len=comtemp_len+mol_len
    if comtemp_len > avr_len:
        comtemp_len=0
        lom=lom+1

    t_sum_len=t_sum_len+mol_len
    i=i+1

print(lom)
print(i)
