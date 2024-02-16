#!/bin/sh
#SBATCH -n 1
#SBATCH --array 1-6
#SBATCH -J testheiretsu
# SBATCH -o testheiretsu.out
# SBATCH -e testheiretsu.err
# COUNT=$(python count.py Compound_119500001_120000000.sdf 500)
COUNT=$1

LAST=$2
echo $SLURM_ARRAY_TASK_ID
echo $COUNT
echo $LAST

start=$(date +%s)
# TIME1=$(cat /proc/uptime | awk '{print $1}')
# echo time1: $TIME1
#python auto_docking.py "PubChem_c_${SLURM_ARRAY_TASK_ID}.sdf" protein.pdb $COUNT "${SLURM_ARRAY_TASK_ID}"
python auto_docking.py "${LAST}" protein.pdb "${COUNT}" "${SLURM_ARRAY_TASK_ID}"
python henkan.py "result${SLURM_ARRAY_TASK_ID}.pdbqt"
end=$(date +%s)

echo "alltime : $(($end-$start)) second"
# TIME2=$(cat /proc/uptime | awk '{print $1}')
# echo time2: $TIME2

# DIFF=$(echo "$TIME2 - $TIME1" | bc)

# echo diff: $DIFF

# python auto_docking.py "PubChem_c_${SLURM_ARRAY_TASK_ID}.sdf" protein.pdb 9 "${SLURM_ARRAY_TASK_ID}"
