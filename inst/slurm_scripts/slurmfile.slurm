#!/bin/sh
#SBATCH --job-name=slurm_sapply
#SBATCH --time=01:10:00
#SBATCH --account=notchpeak-shared-short
#SBATCH --partition=notchpeak-shared-short

module load R/4.2.2
Rscript --vanilla slurmfile.R
