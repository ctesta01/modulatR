#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --array=1-50
#SBATCH --time=12:00:00
#SBATCH --mem-per-cpu=20GB
#SBATCH --job-name=lmtp_simstudy1
#SBATCH --error=sbatch-out/%x_%A_%a.err
#SBATCH --output=sbatch-out/%x_%A_%a.out

module load R

srun Rscript 99_run_everything.R ${SLURM_ARRAY_TASK_ID}
