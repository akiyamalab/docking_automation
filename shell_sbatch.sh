#!/bin/sh
# COUNT='python count.py PubChem_compound_cache_sPcUvdmXvCuLBTQctmR9Nz0ydlJbKAuAcaUQzGq0As1qrT4_records.sdf 9'
COUNT=$(python count.py input_20.sdf 6)

date
P0=$(python convert_sdf_to_pdb.py input_20.sdf "$COUNT" 6)
date
echo $P0
P0=$(python idou.py "$P0")
echo $P0
#sbatch --nodelist=ghost4 slurm.sh $COUNT
#start=$(date +%s)
sbatch  slurm.sh $COUNT "$P0"
#end=$(date +%s)
#echo "all time: $(($end-$start)) second"
