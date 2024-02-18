import sys

'''original=sys.stdout
sys.stdout=open('/dev/null','w')
originale=sys.stderr
sys.stderr=open('/dev/null','w')'''


from openbabel import pybel
import subprocess

#SDF
#protein.pdb

n_cmpds_per_job=int(sys.argv[2])

input_sdf=sys.argv[1]
molecules = pybel.readfile("sdf", input_sdf)
n_cmpds=0

cmd = "prepare_receptor -r protein.pdb -A bonds_hydrogen"
subprocess.run(cmd, shell=True)
    # 複数のPDBファイルに変換
output_pdb_prefix="output"



for idx, mol in enumerate(molecules):
#    output_pdb_file = f"{output_pdb_prefix}_{(idx+(njob-1)*ybun) + 1}.pdb"
    output_pdb_file = f"{output_pdb_prefix}_{idx+ 1}.pdb" 
   #mol.write("pdb", "output_pdb/"+output_pdb_file)#,overwrite=True,opt={"d": None,"h": None,"g": None,"t": None})
    mol.make3D(forcefield="mmff94")
    mol.write("pdb", output_pdb_file,overwrite=True)
    n_cmpds=n_cmpds+1

'''sys.stdout=original
sys.stderr=originale'''

# MEMO: pp? over? (Yanagisawa)
if n_cmpds%(n_cmpds_per_job+1) == 0:
    pp=int(n_cmpds/(n_cmpds_per_job+1))+1
    over=0
else:
    pp=int(n_cmpds/(n_cmpds_per_job+1))+1
    over=n_cmpds-(n_cmpds_per_job*int(n_cmpds/(n_cmpds_per_job+1)))

print(str(pp)+'/'+str(over))
