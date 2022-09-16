#!/bin/bash
# 
#SBATCH -p defq
module load Python
module load matplotlib/3.0.0-foss-2018b-Python-3.6.6

mkdir ../data/${1}saved

python3 imbal_Nr_perc.py -rp $1 -ns $2 -st $3
