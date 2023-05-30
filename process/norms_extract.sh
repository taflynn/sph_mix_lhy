#!/bin/bash
# 
#SBATCH -p defq
module load Python

python3 norms_vs_t.py -rp $1 -sn $2
