#!/bin/bash
# 
#SBATCH -p defq
module load Python

mkdir ../data/${1}saved

python3 iter_norms_vs_omg.py -rp $1 -ns $2
