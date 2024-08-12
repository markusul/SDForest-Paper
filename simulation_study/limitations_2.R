args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[1]))

source('simulation_study/utils.r')

n <- 500
p <- 2
q <- 1
n_test <- 500

N_rep <- 2

seq_ <- seq(1, 1, 3)

performance_measure <- function(n, p, q, n_test, eff, fixEff = FALSE){
  data <- simulate_data_nonlinear(q = q, p = p, n = n+n_test, m = 1, eff = eff, fixEff = fixEff)
  data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
  data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])
  
  fit <- SDForest(Y ~ ., data_train, cp = 0, gpu = F, mc.cores = 10)
  fit2 <- ranger(Y ~ ., data_train, num.trees = 100, importance = 'impurity', mtry = floor(0.5 * ncol(data_train)))
  
  pred <- predict(fit, data_test)
  pred2 <- predict(fit2, data_test)$predictions
  
  mse <- (data$f_X[(n+1):(n+n_test)] - pred)
  mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)
  
  return(list(SDF = mse, ranger = mse2))
}

print('start')
start <- Sys.time()
perf <- lapply(1:N_rep, function(i) lapply(seq_, function(eff) performance_measure(n, p, q, n_test, eff, TRUE)))
save(perf, seq_, file = paste("simulation_study/results/perf_limitations_2/", '1', '.RData', sep=''))
print('done')
print(Sys.time() - start)