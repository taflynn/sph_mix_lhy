#!/bin/bash
# 
#SBATCH -p defq
module load Python

mkdir ../data/${1}saved

python3 decay_fit.py -rp $1 -ns $2 -st $3
