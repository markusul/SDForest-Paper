args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[1]))

source('simulation_study/utils.R')

p <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq_ <- seq(100, 1000, 200)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq_, function(n) performance_measure(n = n, p = p, q = q, n_test = n_test, eff = NULL)))
save(perf, seq_, file = paste("simulation_study/results/perf_n/", args[1], '.RData', sep=''))
print('n done')
print(Sys.time() - start)