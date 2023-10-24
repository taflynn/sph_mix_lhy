#!/bin/bash
# 
#SBATCH -p defq

#SBATCH -J array_centre_saver

#SBATCH -c 1
module load Python

python3 centre_dens_saver.py -rp $1 -jn ${SLURM_ARRAY_TASK_ID} -fn $2
