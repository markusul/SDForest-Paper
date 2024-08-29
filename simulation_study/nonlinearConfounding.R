library(EQL)
library(SDForest)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(tidyr)

mc.cores <- 100

f_hermite <- function(x, beta, js){
    # function to generate f_X
    # x: covariates
    # beta: parameter vector
    # js: relevant covariates

    # number of relevant covariates
    m <- length(js)

    # complexity of f_X
    complexity <- length(beta) / m

    if(is.null(dim(beta))) beta <- matrix(beta)

    # calculate f_X
    do.call(sum, lapply(1:m, function(i) {
        j <- js[i]
        # select beta for covariate j
        res <- lapply(1:complexity, function(k) {
            beta[k, i] * hermite(x[j], k-1)
        })
        Reduce('+', res)
    }))
}

simulate_nonlinear_confounding <- function(q, p, n, m, df){
  H <- matrix(rnorm(n*q), ncol=q)
  
  Betas <- replicate(p, matrix(runif(q*df, -1, 1), nrow = df), simplify = F)
  X <- sapply(Betas, function(beta) {apply(H, 1, function(h) f_hermite(h, beta, 1:q))})
  X <- X + matrix(rnorm(n*p), ncol = p)

  complexity <- 2
  # random parameter for fourier basis
  beta <- runif(m * complexity * 2, -1, 1)
  # random sparse subset of covariates in X
  js <- sample(1:p, m)

  # random coefficient vector delta
  delta <- rnorm(q, 0, 2)

  # generate f_X
  f_X <- apply(X, 1, function(x) f_four(x, beta, js))
  Y <- f_X + H %*% delta + rnorm(n, 0, 0.01)
  return(list(X = X, Y = Y, f_X = f_X, j = js))
}

# Example simulations
set.seed(42)
p <- 30
q <- 1
n <- 50
df <- 4
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
  data <- simulate_nonlinear_confounding(q = q, p = p, n = n+n_test, m = 4, df = 4)
  data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
  data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])
  
  fit <- SDForest(Y ~ ., data_train, mc.cores = mc.cores)
  fit2 <- SDForest(Y ~ ., data_train, Q_type = "no_deconfounding", mc.cores = mc.cores)
  
  pred <- predict(fit, data_test)
  pred2 <- predict(fit2, data_test)
  
  mse <- (data$f_X[(n+1):(n+n_test)] - pred)
  mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)

  return(list(SDF = mse, RF= mse2))
}

perf <- replicate(10, sapply(performance_measure(n, p, q, 500), 
                              function(x)mean(x**2)))
perf <- t(perf)

perf_g <- gather(data.frame(perf), key = 'method', value = 'performance')
save(perf_g, fit, fit2, dep, dep2, sing, df, js, X, Y, f_X,
     file = "simulation_study/results/nonlin_confounding.RData")