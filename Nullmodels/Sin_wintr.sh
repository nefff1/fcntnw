#!/bin/bash

#SBATCH --job-name=sin_wintr

#SBATCH --time=25:00:00
#SBATCH -A node
#SBATCH -p node
#SBATCH --qos=normal
#SBATCH --output=XXX
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4 # 20 is the max
#SBATCH --mem-per-cpu=3000
#SBATCH --array=1-100%50

SCRIPT=wintr_networkanalyses.R

srun singularity exec ubuntu_R3.5.2_ver2 Rscript --vanilla $SCRIPT $SLURM_ARRAY_TASK_ID "raw" "T" "F" "wintr_data.Rdata" 2999




