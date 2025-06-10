module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
export OMP_NUM_THREADS=1

for i in {1..20}
do 
   sbatch --time=10:00:00 --job-name="n $i" --mem-per-cpu=7GB --output=outfiles/n_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_n.R $i"
   sbatch --time=10:00:00 --job-name="p $i" --mem-per-cpu=4GB --output=outfiles/p_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_p.R $i"
   sbatch --time=10:00:00 --job-name="q $i" --mem-per-cpu=4GB --output=outfiles/q_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_q.R $i"
   sbatch --time=10:00:00 --job-name="eff $i" --mem-per-cpu=4GB --output=outfiles/eff_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/dense_assumption.R $i"
done