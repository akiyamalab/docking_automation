import time
import subprocess


from openbabel import openbabel
from openbabel import pybel
import sys

import os


procount=int(sys.argv[5])


ybun=int(sys.argv[3])+1
njob=int(sys.argv[4])

fg=open('joblist_'+str(njob)+'.txt','r')
fglist=fg.readlines()


cmd="touch taioulist"+str(njob)+".txt"
subprocess.run(cmd, shell=True)
taioulist_file=open("taioulist"+str(njob)+".txt",'w')

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

#sys.argv[1]=sdf
#sys.argv[2]=pdb
#sys.argv[3]=pdbconf
"""def convert_sdf_to_pdb(input_sdf, output_pdb_prefix):
    molecules = pybel.readfile("sdf", input_sdf)
    i=0
    for idx, mol in enumerate(molecules):
        output_pdb_file = f"{output_pdb_prefix}_{(idx+(njob-1)*ybun) + 1}.pdb"
        #mol.write("pdb", "output_pdb/"+output_pdb_file)#,overwrite=True,opt={"d": None,"h": None,"g": None,"t": None})
        #print(i)
        mol.make3D(forcefield="mmff94")
        mol.write("pdb", output_pdb_file,overwrite=True)
        i=i+1
    return i"""

if __name__ == "__main__":
    time1=time.time()
#    input_sdf_file = sys.argv[1]
#    n=int(sys.argv[1])
    output_pdb_prefix = "output"

    print("aaaa="+sys.argv[1]+'aaaaaaaa')

    #last=int(sys.argv[1].split('\n')[2].split('/')[1])
    #point=int(sys.argv[1].split('\n')[2].split('/')[0])
    last=int(sys.argv[1].split('/')[1])
    point=int(sys.argv[1].split('/')[0])

    if njob < point:
        n=ybun
    elif njob == point:
        n=last
    else:
        n=0

    n=len(fglist)

    print('xxxxxxxxxxx')

    #n=convert_sdf_to_pdb(input_sdf_file, output_pdb_prefix)
    #cmd="rm "+input_sdf_file
    #subprocess.run(cmd, shell=True)

    #cmd = "conda deactivate"
    #subprocess.run(cmd, shell=True)
    #command = f'''
    #load protein.pdb
    #remove solvent
    #save protein.pdb
    #'''
    #subprocess.run([pymol_path, '-c', '-q', '-d', command], check=True)
    #for z in range(n):
    #    remove_water_from_pdb('output_'+str(z+1)+'.pdb', 'output_'+str(z+1)+'.pdb')"""


#    cmd = "prepare_receptor -r "+sys.argv[2]+" -A bonds_hydrogen"
#    subprocess.run(cmd, shell=True)
    time2=time.time()
    #print("prepare time="+str(time2-time1))
    time1=time.time()


    for pnow in procount:
        cmd="> result"+str(pnow)+"/result"+str(njob)+".pdbqt"
        subprocess.run(cmd, shell=True)
        cmd = "mkdir result"+str(pnow)
        subprocess.run(cmd, shell=True)

        input_taiou_file=open('data_taiou.txt','r')
        taiou_data_list=input_taiou_file.readlines()
        input_taiou_file.close()
        njob_taiou_list=open('njob_taiou_'+str(njob)+'.txt','w')
        allpre_taiou=open('data_taiou.txt','r')
        allpre_taiou_list=allpre_taiou.readlines()
        allpre_taiou.close()

        for e in range(n):
            filename=fglist[e].replace('\n','')
            if filename == '':
                break
            input_file=open(filename,'r')
            lines=input_file.readlines()
            for line in lines:
                if line[0:12]=='REMARK  Name':
                    title=line.split()[3]
            input_file.close()
            outfilename='result'+str(pnow)+'/tto'+str(njob)+'_'+str(e+1)+'.txt'
            cmd='touch '+outfilename
            subprocess.run(cmd, shell=True)
            #filename='output_'+str((e+(njob-1)*ybun)+1)+'.pdbqt'
            if not os.path.exists(filename):
                continue
            weight = get_molecule_weight(filename)
            if weight >= 700:
                for tai in allpre_taiou_list:
                    if tai == '':
                        break
                    if tai.split(',')[0] == title:
                        njob_taiou_list.write(tai.replace('\n','')+',C'+'\n')
                continue
            #print((e+(njob-1)*ybun)+1)
            time1=time.time()
            cmd="~/autodock_vina_1_1_2_linux_x86/bin/vina --seed 42 --cpu 1 --num_modes 1 --receptor input_protein/protein"+str(pnow)+".pdbqt --ligand "+filename+" --config input_conf/autodock"+str(pnow)+".conf --out result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".pdbqt --log result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".log > "+outfilename
            subprocess.run(cmd, shell=True)
            time2=time.time()
            dockingtime=time2-time1

            for tai in allpre_taiou_list:
                if tai == '':
                    break
                if tai.split(',')[0] == title:
                    njob_taiou_list.write(tai.replace('\n','')+','+str(dockingtime)+'\n')
                    break

            #cmd="rm output_"+str((e+(njob-1)*ybun)+1)+".pdbqt"
            #subprocess.run(cmd, shell=True)
            cmd="touch result"+str(pnow)+"/result"+str(njob)+".pdbqt"
            subprocess.run(cmd, shell=True)
            cmd="chmod 755 result"+str(pnow)+"/result"+str(njob)+".pdbqt"
            subprocess.run(cmd, shell=True)
            cmd="chmod 755 result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".pdbqt"
            subprocess.run(cmd, shell=True)
            cmd="cat result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".pdbqt >> result/result"+str(njob)+".pdbqt"         #追記する
            subprocess.run(cmd, shell=True)
            cmd="rm result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".pdbqt"            #削除する
            subprocess.run(cmd, shell=True)
            cmd="rm result"+str(pnow)+"/multi_autodock"+str(njob)+'_'+str(e+1)+".log"            #削除する
            subprocess.run(cmd, shell=True)

            #taioulist_file.write(str((e+(njob-1)*ybun)+1)+' ')
        time2=time.time()

        taioulist_file.close()
        njob_taiou_list.close()
        #print("docking time="+str(time2-time1))

