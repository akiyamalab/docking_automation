#!/bin/sh
#SBATCH -n 1
#SBATCH --array 1-10
#SBATCH -J testheiretsu
# SBATCH -o testheiretsu.out
# SBATCH -e testheiretsu.err
# COUNT='python count.py PubChem_compound_cache_sPcUvdmXvCuLBTQctmR9Nz0ydlJbKAuAcaUQzGq0As1qrT4_records.sdf 9'
echo $SLURM_ARRAY_TASK_ID
# python auto_docking.py "PubChem_c_${SLURM_ARRAY_TASK_ID}.sdf" protein.pdb "${COUNT}" "${SLURM_ARRAY_TASK_ID}"
python auto_docking.py "PubChem_c_${SLURM_ARRAY_TASK_ID}.sdf" protein.pdb 10 "${SLURM_ARRAY_TASK_ID}"
