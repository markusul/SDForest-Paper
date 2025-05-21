module load stack/2024-06  gcc/12.2.0
module load r/4.4.0
export OMP_NUM_THREADS=1

for i in {1..12}
do 
   sbatch --time=128:00:00 --job-name="n $i" --mem-per-cpu=2GB --output=outfiles/n_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_n.r $i"
   sbatch --time=48:00:00 --job-name="p $i" --mem-per-cpu=2GB --output=outfiles/p_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_p.r $i"
   sbatch --time=128:00:00 --job-name="q $i" --mem-per-cpu=2GB --output=outfiles/q_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/perf_q.r $i"
   sbatch --time=128:00:00 --job-name="eff $i" --mem-per-cpu=2GB --output=outfiles/eff_$i.out --cpus-per-task=100 --wrap "Rscript --vanilla simulation_study/dense_assumption.r $i"
done