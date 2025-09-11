#!/bin/bash
#SBATCH --time=10:00:00
#SBATCH --job-name="PlainSDF"
#SBATCH --mem-per-cpu=2GB
#SBATCH --output=outfiles/PlainSDF.out
#SBATCH --cpus-per-task=100

module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
module load python/3.11.6
export OMP_NUM_THREADS=1

Rscript --vanilla cBench/R/PlainSDF.R