import time
import subprocess
#cmd = "export PYTHONPATH=/Users/kaito/mgltools_1.5.7_MacOS-X/MGLToolsPckgs/:$PYTHONPATH"
#subprocess.run(cmd, shell=True)
#cmd = "conda activate tmp_20231111_2"
#subprocess.run(cmd, shell=True)

from openbabel import openbabel
from openbabel import pybel
import sys
#sys.path.append("/usr/local/bin/pymol")
#from pymol import cmd
import os

print(sys.argv[3])
ybun=int(sys.argv[3])+1
njob=int(sys.argv[4])




if __name__ == "__main__":
    time1=time.time()
#    input_sdf_file = sys.argv[1]
#    n=int(sys.argv[1])
    output_pdb_prefix = "ligand/output"

    print("aaaa="+sys.argv[1]+'aaaaaaaa')

#    last=int(sys.argv[1].split('\n')[2].split('/')[1])
#    point=int(sys.argv[1].split('\n')[2].split('/')[0])
    last=int(sys.argv[1].split('/')[1])
    point=int(sys.argv[1].split('/')[0])

    if njob < point:
        n=ybun
    elif njob == point:
        n=last
    else:
        n=0



    if njob <= point:
        cmd='mkdir Pub'+str(njob)
        subprocess.run(cmd, shell=True)
        cmd='mk_prepare_ligand_2.py -i PubChem_c_'+str(njob)+'.sdf --multimol_outdir Pub'+str(njob)+' --multimol_prefix output_'
        subprocess.run(cmd, shell=True)

        for hum in range(n):
            cmd='cp Pub'+str(njob)+'/output_-'+str(hum+1)+'.pdbqt output_'+str(ybun*(njob-1)+hum+1)+'.pdbqt'
            subprocess.run(cmd, shell=True)
            cmd='rm Pub'+str(njob)+'/output_-'+str(hum+1)+'.pdbqt'
            subprocess.run(cmd, shell=True)

    print('xxxxxxxxxxx')


"""    for h in range(n):
        print(h)
        #print((h+(njob-1)*ybun)+1)
        #filename="/mnt/fs/k_takahashi/test_20240422/ligand/output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        #filename="output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        filename='output_'+str((h+(njob-1)*ybun)+1)+'.pdb'
        if os.path.exists(filename):
           print("sonzaisuru")
        cmd="cd ligand"
        #subprocess.run(cmd, shell=True)

        cmd = "prepare_ligand -l "+filename+" -o "+filename+"qt -A bonds_hydrogen"
        try:
           subprocess.run(cmd, shell=True)
        except subprocess.TimeoutExpired:
          # cmd="rm ligand/output_"+str((h+(njob-1)*ybun)+1)+".pdb"
         #  cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
           cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
           subprocess.run(cmd, shell=True)
           continue
        #cmd="rm ligand/output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        #cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        cmd="rm output_"+str((h+(njob-1)*ybun)+1)+".pdb"
        subprocess.run(cmd, shell=True)
        cmd="cd .."
        #subprocess.run(cmd, shell=True)
"""

