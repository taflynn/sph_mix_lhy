#!/bin/bash
# 
#SBATCH -t 00:05:00
#SBATCH -p defq
module load Python
module load matplotlib/3.2.1-foss-2020a-Python-3.8.2
python3 running.py
