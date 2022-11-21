#!/bin/bash
# 
#SBATCH -p defq
module load Python
module load matplotlib/3.0.0-foss-2018b-Python-3.6.6

mkdir ../data/${1}output

python3 mu_width_acs_scan.py -rp $1 -ns $2
