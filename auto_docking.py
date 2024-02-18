print('fffffffff')
import time
import subprocess
#cmd = "export PYTHONPATH=/Users/kaito/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/:$PYTHONPATH"
#subprocess.run(cmd, shell=True)
#cmd = "conda activate tmp_20231111_2"
#subprocess.run(cmd, shell=True)

import sys
#sys.path.append("/usr/local/bin/pymol")
#from pymol import cmd
import os

import configparser
config = configparser.ConfigParser()
config.read('settings.ini')
VINA             = config["EXECUTABLE"]["VINA"]
PREPARE_LIGAND   = config["EXECUTABLE"]["PREPARE_LIGAND"]

ybun=int(sys.argv[3])+1
njob=int(sys.argv[4]) # n-th job (not number of jobs)


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
    
    last=int(sys.argv[1].split('\n')[2].split('/')[1])
    point=int(sys.argv[1].split('\n')[2].split('/')[0])

    if njob < point:
        n=ybun
    elif njob == point:
        n=last
    else:
        n=0  
  
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

    for h in range(n):
        print(h)
        #print((h+(njob-1)*ybun)+1)
        filename='output_'+str((h+(njob-1)*ybun)+1)+'.pdb'
        cmd = f"{PREPARE_LIGAND} -l {filename} -A bonds_hydrogen"
        try:
           subprocess.run(cmd, shell=True)
        except subprocess.TimeoutExpired:
           cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
           subprocess.run(cmd, shell=True)
           continue
        cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        subprocess.run(cmd, shell=True)
#    cmd = "prepare_receptor -r "+sys.argv[2]+" -A bonds_hydrogen"
#    subprocess.run(cmd, shell=True)
    time2=time.time()
    print("prepare time="+str(time2-time1))
    time1=time.time()
    cmd="> result/result"+str(njob)+".pdbqt"
    subprocess.run(cmd, shell=True)
    cmd = "mkdir result"
    subprocess.run(cmd, shell=True)
    for e in range(n):
        outfilename='result/tto'+str((e+(njob-1)*ybun)+1)+'.txt'
        cmd='touch '+filename
        subprocess.run(cmd, shell=True)
        filename='output_'+str((e+(njob-1)*ybun)+1)+'.pdbqt'
        if not os.path.exists(filename):
            continue
        #print((e+(njob-1)*ybun)+1)
        cmd=f"{VINA} --seed 42 --cpu 1 --num_modes 1 --receptor protein.pdbqt --ligand "+filename+" --config autodock.conf --out result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".pdbqt --log result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".log > "+outfilename
        subprocess.run(cmd, shell=True)
        cmd="rm output_"+str((e+(njob-1)*ybun)+1)+".pdbqt"
        subprocess.run(cmd, shell=True)
        cmd="touch result/result"+str(njob)+".pdbqt"
        subprocess.run(cmd, shell=True)
        cmd="chmod 755 result/result"+str(njob)+".pdbqt"
        subprocess.run(cmd, shell=True)
        cmd="chmod 755 result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".pdbqt"
        subprocess.run(cmd, shell=True)
        cmd="cat result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".pdbqt >> result/result"+str(njob)+".pdbqt"         #追記する
        subprocess.run(cmd, shell=True)
        cmd="rm result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".pdbqt"            #削除する
        subprocess.run(cmd, shell=True)
        cmd="rm result/multi_autodock"+str((e+(njob-1)*ybun)+1)+".log"            #削除する
        subprocess.run(cmd, shell=True)
    time2=time.time()

    print("docking time="+str(time2-time1))
