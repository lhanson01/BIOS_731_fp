#!/bin/bash
#SBATCH --array=1-4
#SBATCH --job-name=liam_731_bgl_real
#SBATCH --partition=wrobel
#SBATCH --output=bgl_real.out
#SBATCH --error=bgl_real.err

module purge
module load R

JOBID=$SLURM_ARRAY_TASK_ID
Rscript source/03_fit_real_data.R $JOBID