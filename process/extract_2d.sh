#!/bin/bash
# 
#SBATCH -p defq
module load Python

mkdir -p ../data/${1}saved

python3 proc_2d_omg_imb.py -rp $1 -on $2 -in $3
