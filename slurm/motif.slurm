#!/bin/bash
#SBATCH --time=350:00:00
#SBATCH --job-name="motif"
#SBATCH --mem-per-cpu=100GB
#SBATCH --output=outfiles/motif.out
#SBATCH --cpus-per-task=5

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla motif/motif_fit.R