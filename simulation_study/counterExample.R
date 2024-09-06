library(SDForest)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(tidyr)
library(parallel)

mc.cores <- 100


simulate_data_counter <- function(q, p, n){
  # random confounding covariates H
  H <- matrix(rnorm(n * q, 0, 1), nrow = n)
  
  beta <- runif(q, -1, 1)
  
  # random correlation matrix cov(X, H)
  Gamma <- matrix(rnorm(q * p, 0, 1), nrow = q)
  
  # random coefficient vector delta
  delta <- rnorm(q, 0, 1)
  
  # random error term
  E <- matrix(rnorm(n * p, 0, 1), nrow = n)
  
  if(q == 0){
    X <- E
  }else{
    X <- H %*% Gamma + E
  }
  
  f_X <- X %*% t(Gamma ^ (-1)) %*% beta
  
  # generate Y
  Y <- f_X + H %*% delta + rnorm(n, 0, 0.1)
  
  #return data
  list(X = X, Y = Y, f_X = f_X, beta = beta, H = H)
}


# Example simulations
set.seed(42)
p <- 300
q <- 1
n <- 500

# Quantitative performance comparison
performance_measure <- function(n, p, q, n_test){
  data <- simulate_data_counter(q = q, p = p, n = n+n_test)
  data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
  data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])
  
  fit <- SDForest(Y ~ ., data_train, mc.cores = 1)
  pred <- predict(fit, data_test)
  
  fit <- SDForest(Y ~ ., data_train, Q_type = "no_deconfounding", mc.cores = 1)
  pred2 <- predict(fit, data_test)
  
  mse <- (data$f_X[(n+1):(n+n_test)] - pred)
  mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)
  
  return(list(SDF = mse, RF= mse2))
}

perf <- mclapply(1:100, function(i) sapply(performance_measure(n, p, q, 500), 
                              function(x)mean(x**2)), mc.cores = mc.cores)
perf <- do.call(rbind, perf)

perf_g <- gather(data.frame(perf), key = 'method', value = 'performance')
save(perf_g,
     file = "simulation_study/results/counterExample.RData")
