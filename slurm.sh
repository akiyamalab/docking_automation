#!/bin/sh
#SBATCH -n 1
#SBATCH --array 1-352
#SBATCH -J testheiretsu

COUNT=$1

LAST=$2
echo $SLURM_ARRAY_TASK_ID
echo $COUNT
echo $LAST

start=$(date +%s)

python auto_docking.py "${LAST}" protein.pdb "${COUNT}" "${SLURM_ARRAY_TASK_ID}"
#python henkan.py "result${SLURM_ARRAY_TASK_ID}.pdbqt" "${SLURM_ARRAY_TASK_ID}"
end=$(date +%s)

echo "alltime : $(($end-$start)) second"

