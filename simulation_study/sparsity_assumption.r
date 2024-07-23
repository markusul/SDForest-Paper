args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[1]))

source('simulation_study/utils.r')

n <- 500
p <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq <- seq(0, 499, 50)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq, function(eff) performance_measure(n, p, q, n_test, eff)))
save(perf, seq, file = paste("simulation_study/results/perf_eff/", args[1], '.RData', sep=''))
print('done')
print(Sys.time() - start)