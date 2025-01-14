library(SDModels)
library(parallel)
library(tidyr)

mc.cores <- 100

performance_measure <- function(n, p, q, n_test){
  data <- simulate_data_step(q = q, p = p, n = n+n_test, m = 10)
  data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
  data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])
  
  fit <- SDTree(Y ~ ., data_train)
  pred <- predict(fit, data_test)
  
  fit <- SDTree(Y ~ ., data_train, fast = FALSE)
  pred2 <- predict(fit, data_test)
  
  mse <- (data$f_X[(n+1):(n+n_test)] - pred)
  mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)
  
  return(list(SDT1 = mse, SDT2= mse2))
}


perf <- mclapply(1:100, function(i) sapply(performance_measure(1000, 500, 20, 500), 
                                           function(x)mean(x**2)), mc.cores = mc.cores)
perf <- do.call(rbind, perf)
perf_g <- gather(data.frame(perf), key = 'method', value = 'performance')

perf_g
save(perf_g,
     file = "simulation_study/results/fast.RData")