#!/bin/bash
#SBATCH --time=350:00:00
#SBATCH --job-name="cBench"
#SBATCH --mem-per-cpu=30GB
#SBATCH --output=outfiles/cBench_sub.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
module load python/3.11.6

export OMP_NUM_THREADS=1

Rscript --vanilla cBench/R/fit.R