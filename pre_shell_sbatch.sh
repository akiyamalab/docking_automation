#!/bin/sh
# COUNT='python count.py PubChem_compound_cache_sPcUvdmXvCuLBTQctmR9Nz0ydlJbKAuAcaUQzGq0As1qrT4_records.sdf 9'
mkdir ligand
chmod 755 ligand

rm output*.pdbqt

PROCOUNT=$(python receptor_N_count.py)
for ((N=1; N<=PROCOUNT; N++))
do
  prepare_receptor -r input_protein/protein${N}.pdb -o input_protein/protein${N}.pdbqt -A bonds_hydrogen
done

#prepare_receptor -r protein.pdb -A bonds_hydrogen

COUNT=$(python count.py input_100000.sdf 352)

date
# P0=$(python pre_convert_sdf_to_pdb.py input_100000.sdf "$COUNT" 352)
#python pre_convert_sdf_to_pdb.py input_100000.sdf "$COUNT" 352

python pre_devide.py input_100000.sdf 352

P0=$(python dock_P0.py input_100000.sdf "$COUNT" 352)
date
echo $P0
P0=$(python idou.py "$P0")
echo $P0


sbatch  pre_slurm.sh $COUNT "$P0"
