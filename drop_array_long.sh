#!/bin/bash
#
# set up a job array with tasks numbered 1 to n.
#
#SBATCH -p long
#SBATCH -t 14-00:00:00
# give the array a single job name
#SBATCH -J array_sph_drop
#
# request one core per task
#SBATCH -c 1
#
module load Python
module load matplotlib/3.0.0-foss-2018b-Python-3.6.6
#
# use the task ID to locate the input file for the task.
python3 running.py -wp $1${SLURM_ARRAY_TASK_ID} -rp ./config_dens_$2${SLURM_ARRAY_TASK_ID}.json
