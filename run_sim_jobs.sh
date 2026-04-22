#!/bin/bash
#SBATCH --array=1-400%100
#SBATCH --job-name=liam_731_bgl_sim
#SBATCH --partition=wrobel
#SBATCH --output=bgl_sim.out
#SBATCH --error=bgl_sim.err

module purge
module load R

JOBID=$SLURM_ARRAY_TASK_ID
Rscript source/run_sim_job.R $JOBID