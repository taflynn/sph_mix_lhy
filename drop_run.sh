#!/bin/bash
# 
#SBATCH -t 00:05:00
#SBATCH -p defq
#SBATCH --ntasks = 1
#SBATCH --cpus-per-task = 1
#SBATCH --ntasks-per-core = 1
module load Python
python3 running.py
