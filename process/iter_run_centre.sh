#!/bin/bash
# 
#SBATCH -p defq
module load Python

python3 iter_centre_dens_saver.py -rp $1 -nj $2
