#!/bin/bash
# 
#SBATCH -p defq
module load Python

python3 centre_dens_saver.py -rp $1 -jn $2 -fn $3
