#!/bin/bash
# 
#SBATCH -p long
#SBATCH -t 14-00:00:00
module load Python
module load matplotlib/3.0.0-foss-2018b-Python-3.6.6
python3 running.py -wp $1 -rp $2
