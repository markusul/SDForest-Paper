#!/bin/bash
#SBATCH --time=50:00:00
#SBATCH --job-name="fast"
#SBATCH --mem-per-cpu=1GB
#SBATCH --output=outfiles/fast.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0

export OMP_NUM_THREADS=1

Rscript --vanilla simulation_study/fast.R