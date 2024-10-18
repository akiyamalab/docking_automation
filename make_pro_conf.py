#引数はoriginalフォルダ（と同等のフォルダ）のパス

import numpy as np
import sys
import subprocess
import os

def parse_pdb(file):
    atoms = []
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                # X, Y, Z 座標を取得
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append([x, y, z])
    return np.array(atoms)

def calculate_center_of_mass(atoms):
    # 各軸の座標の平均を計算
    return np.mean(atoms, axis=0)

original_path=sys.argv[1]

prolist=os.listdir(original_path)
p=1
outputfile=open("pro_number_log.txt",'w')

cmd="mkdir input_protein"
subprocess.run(cmd, shell=True)
cmd="chmod 777 input_protein"
subprocess.run(cmd, shell=True)
cmd="mkdir input_conf"
subprocess.run(cmd, shell=True)
cmd="chmod 777 input_conf"
subprocess.run(cmd, shell=True)

for pname in prolist:
    outputfile.write(str(p)+' '+pname+'\n')
    cmd="cp "+original_path+"/"+pname+"/"+"receptor.pdb"+" "+"input_protein/protein"+str(p)+".pdb"
    subprocess.run(cmd, shell=True)
    cmd="chmod 755 input_protein/protein"+str(p)+".pdb"
    subprocess.run(cmd, shell=True)
    cmd="./eBox.pl "+original_path+"/"+pname+"/crystal_ligand.mol2"
    #result=subprocess.run(cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True)
    result=subprocess.run(cmd, shell=True,stdout=subprocess.PIPE,text=True)
    outresult=result.stdout
    print(pname)
    #print(outresult)
    if outresult == '':
        p=p+1
        continue
    size=outresult.split('\n')[1]


    recep_center=calculate_center_of_mass(parse_pdb(original_path+"/"+pname+"/"+"receptor.pdb"))
    cmd="touch input_conf/autodock"+str(p)+".conf"
    subprocess.run(cmd, shell=True)
    conffile=open("input_conf/autodock"+str(p)+".conf",'w')
    conffile.write("center_x = "+str(recep_center[0])+'\n')
    conffile.write("center_y = "+str(recep_center[1])+'\n')
    conffile.write("center_z = "+str(recep_center[2])+'\n')
    conffile.write("size_x = "+size.split(' ')[0]+'\n')
    conffile.write("size_y = "+size.split(' ')[1]+'\n')
    conffile.write("size_z = "+size.split(' ')[2]+'\n')

    conffile.close()
    p=p+1

outputfile.close()

