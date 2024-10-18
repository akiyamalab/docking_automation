#!/bin/sh
COUNT=$(python count.py input_100000.sdf 352)

python equal_multi.py
PROCOUNT=$(python receptor_N_count.py)

date
P0=$(python dock_P0.py input_100000.sdf "$COUNT" 352)
date
echo $P0
P0=$(python idou.py "$P0")
echo $P0

sbatch  slurm.sh $COUNT "$P0" $PROCOUNT
