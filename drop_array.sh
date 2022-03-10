#!/bin/bash
#
# set up a job array with tasks numbered 1 to 10.
#SBATCH --array=1-10
#
# give the array a single job name
#SBATCH -J array_sph_drop
#
# Standard error messages are saved in this file
#SBATCH -e drop_array_%A_%a.err
#
# request one core per task
#SBATCH -c 1
#
module load Python
module load matplotlib/3.0.0-foss-2018b-Python-3.6.6
#
# use the task ID to locate the input file for the task.
python3 running.py -wp ./data/drop_array_%A_%a -rp ./${SLURM_ARRAY_TASK_ID}.json
