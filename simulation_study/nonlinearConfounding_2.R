library(SDModels)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(tidyr)
library(parallel)

mc.cores <- 100


simulate_nonlinear_confounding <- function(q, p, n, m, df){
  H <- matrix(rnorm(n*q), ncol=q)
  complexity_H <- df
  
  X <- sapply(1:p, function(j){
    beta_H <- runif(q * complexity_H * 2, -1, 1)
    apply(H, 1, function(x) f_four(x, beta_H, 1:q))
  })
  X <- X + matrix(rnorm(n*p, 0, 1), ncol = p)

  complexity <- 2
  # random parameter for fourier basis
  beta <- runif(m * complexity * 2, -1, 1)
  # random sparse subset of covariates in X
  js <- sample(1:p, m)

  # random coefficient vector delta
  delta <- runif(q * complexity_H * 2, -2, 2)
  H_js <- 1:q

  # generate f_X
  f_X <- apply(X, 1, function(x) f_four(x, beta, js))
  f_H <- apply(H, 1, function(x) f_four(x, delta, H_js))
  Y <- f_X + f_H + rnorm(n, 0, 0.01)
  return(list(X = X, Y = Y, f_X = f_X, j = js, H = H))
}

# Example simulations
set.seed(42)
p <- 300
q <- 1
n <- 500
df <- 12
m <- 1

dat <- simulate_nonlinear_confounding(q = q, p = p, n = n, m = m, df = df)
X <- dat$X
Y <- dat$Y
f_X <- dat$f_X
js <- dat$j

# spiking of df singular values
Q <- get_Q(X, 'trim')
d <- svd(X)$d
sing <- data.frame(d, i = 1:p)

# transformed response and f(X)
df <- data.frame(f = f_X, Y, QY = Q %*% Y, Qf = Q %*% f_X)

# Estimation of models
fit <- SDForest(x = X, y = Y, mc.cores = mc.cores)
dep <- partDependence(fit, js[1], mc.cores = mc.cores)
fit2 <- SDForest(x = X, y = Y, Q_type = 'no_deconfounding', mc.cores = mc.cores)
dep2 <- partDependence(fit2, js[1], mc.cores = mc.cores)

# Quantitative performance comparison

performance_measure <- function(n, p, q, n_test){
  data <- simulate_nonlinear_confounding(q = q, p = p, n = n+n_test, m = 4, df = 12)
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
save(perf_g, fit, fit2, dep, dep2, sing, df, js, X, Y, f_X,
     file = "simulation_study/results/nonlin_confounding_2.RData")