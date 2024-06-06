#!/bin/bash
# 
#SBATCH -p defq
module load Python

mkdir ../data/${1}saved

python3 drop_width_saver.py -rp $1 -jn $2
