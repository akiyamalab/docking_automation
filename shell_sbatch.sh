#!/bin/sh
# N_CMPDS_PER_FILE='python count.py PubChem_compound_cache_sPcUvdmXvCuLBTQctmR9Nz0ydlJbKAuAcaUQzGq0As1qrT4_records.sdf 9'
N_CMPDS_PER_FILE=$(python count.py input_20.sdf 6)

date
P0=$(python convert_sdf_to_pdb.py input_20.sdf "$N_CMPDS_PER_FILE" 6)
# Example: P0= "adding gasteiger charges to peptide\n6/12"
date
echo $P0
P0=$(python idou.py "$P0")
# Example: P0= "adding gasteiger charges to peptide\n6/12\n6/12"
echo $P0
#sbatch --nodelist=ghost4 slurm.sh $N_CMPDS_PER_FILE
#start=$(date +%s)
sbatch  slurm.sh $N_CMPDS_PER_FILE "$P0"
#end=$(date +%s)
#echo "all time: $(($end-$start)) second"
