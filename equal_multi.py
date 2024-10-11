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

sum_len=0.0
num=0
original_list=[]

for filename in pdbqt_files:
    input_file=open(filename,'r')
    lines=input_file.readlines()
    for line in lines:
        if line[0:12]=='REMARK  Name':
            num=num+1
            print(line.split())
            title=line.split()[3]
            mol_len=int(get_molecule_length(filename))
            f_len=0.5183*float(mol_len)*float(mol_len)-17.524*float(mol_len)+172.34
            sum_len=sum_len+f_len
            mol_tor=int(count_active_torsions(filename))
            weight=get_molecule_weight(filename)

            output_taiou_file.write(title+','+str(mol_len)+','+str(mol_tor)+','+str(weight)+'\n')

            original_list.append([filename,f_len])



            if float(weight) > 700:
                output_over_file.write(title+'\n')

output_over_file.close()
output_taiou_file.close()

avr_len=sum_len/352

#print('sum_len='+str(sum_len))
#print('avr_len='+str(avr_len))

output_job_file=[]
for j in range(352):
    output_job_file.append(open('joblist_'+str(j+1)+'.txt','w'))

#sublists=[[][],0.0] for _ in range(352)]
sublists=[]
for i in range(352):
    sublists.append([[],0.0])

sorted_list=sorted(original_list,key=lambda x: x[1],reverse=True)

'''for item in sorted_list:
    id_=item[0]
    number=item[1]
    # 合計がNに最も近いサブリストを探す
    min_sum_index = min(range(352), key=lambda i: abs(sublists[i][1] + number - avr_len))
    # サブリストにIDを追加し、数値の合計を更新
    sublists[min_sum_index][0].append(id_)
    #sublists[min_sum_index] = [sublists[min_sum_index][0], sublists[min_sum_index][1] + number]
    sublists[min_sum_index][1] = sublists[min_sum_index][1] + number'''

i=0
for yylis in sorted_list:
    sorted_sublist=sorted(sublists,key=lambda x: x[1])
    j=0
    for sub_y in sorted_sublist:
        if yylis[1]+sub_y[1] < avr_len:
            sub_y[0].append(yylis[0])
            sub_y[1]=sub_y[1]+yylis[1]
            j=1
            break
    if j == 0:
        sorted_sublist[0][0].append(yylis[0])
        sorted_sublist[1][1]=sorted_sublist[1][1]+yylis[1]




print('Average= '+str(avr_len))

# サブリストの合計を出力して確認
i=0
for sub in sorted_sublist:
    for kp in sub[0]:
        output_job_file[i].write(kp+'\n')
    print('Sublist'+str(i+1)+': '+str(sub[1]))
    i=i+1

# オプション: 各サブリストの長さを確認
#    print(f"Sublist {i+1} length: {len(id_list)}")


'''
i=0
comtemp_len=0
t_sum_len=0
lom=0
for filename in pdbqt_files:
    output_job_file[lom].write(filename+'\n')
    mol_len=int(get_molecule_length(filename))
    f_len=0.5183*float(mol_len)*float(mol_len)-17.524*float(mol_len)+172.34
    comtemp_len=comtemp_len+f_len
    if comtemp_len > avr_len:
        comtemp_len=0
        lom=lom+1

    t_sum_len=t_sum_len+mol_len
    i=i+1
print(lom)'''
