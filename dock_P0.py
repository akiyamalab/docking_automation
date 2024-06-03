import sys


from openbabel import openbabel
from openbabel import pybel
import subprocess
import time


c=int(sys.argv[2])
b=int(sys.argv[3])


#pointを超えたらsys.exit()


input_sdf=sys.argv[1]
molecules = pybel.readfile("sdf", input_sdf)
i=0

#cmd = "prepare_receptor -r protein.pdb -A bonds_hydrogen"
#subprocess.run(cmd, shell=True)
    # 複数のPDBファイルに変換
output_pdb_prefix="output"
#output_pdb_prefix="ligand/output"



for idx, mol in enumerate(molecules):
    #output_pdb_file = f"{output_pdb_prefix}_{idx+ 1}.pdb"
    #mol.make3D(forcefield="mmff94")
    #mol.write("pdb", output_pdb_file,overwrite=True)
    i=i+1



if i%(c+1) == 0:
    pp=int(i/(c+1))+1
    over=0
else:
    pp=int(i/(c+1))+1
    over=i-((c+1)*int(i/(c+1)))

print(str(pp)+'/'+str(over))
