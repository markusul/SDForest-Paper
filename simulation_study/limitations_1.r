args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[1]))

source('simulation_study/utils.r')

n <- 500
p <- 500
q <- 20
n_test <- 500

N_rep <- 10

seq_ <- seq(4, 100, 10)

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq_, function(eff) performance_measure(n, p, q, n_test, eff, TRUE)))
save(perf, seq_, file = paste("simulation_study/results/perf_limitations_1/", args[1], '.RData', sep=''))
print('done')
print(Sys.time() - start)