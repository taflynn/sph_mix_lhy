#!/bin/bash
# 
#SBATCH -p defq
module load Python

mkdir ../data/${1}saved_radcut${3}

python3 iter_norms_vs_omg.py -rp $1 -ns $2 -rc $3
