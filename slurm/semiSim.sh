module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
module load python/3.11.6
export OMP_NUM_THREADS=1

for j in {1..20}
do
   for i in {1..11}
   do 
      sbatch --time=128:00:00 --job-name="semi $i $j" --mem-per-cpu=2GB --output=outfiles/semi_$i-$j.out --cpus-per-task=100 --wrap "Rscript --vanilla cBench/R/semiSim.R $i $j"
   done
done