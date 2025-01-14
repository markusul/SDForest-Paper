args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[1]))

source('simulation_study/utils.R')

n <- 500
p <- 500

n_test <- 500

N_rep <- 10

seq_ <- seq(0, 100, 20)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq_, function(q) performance_measure(n = n, p = p, q = q, n_test = n_test, eff = NULL)))
save(perf, seq_, file = paste("simulation_study/results/perf_q/", args[1], '.RData', sep=''))
print('q done')
print(Sys.time() - start)
