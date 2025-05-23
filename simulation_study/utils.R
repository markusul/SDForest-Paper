library(SDModels)
library(ranger)
library(parallel)

mc.cores <- 100

performance_measure <- function(n, p, q, n_test, eff, fixEff = FALSE){
    data <- simulate_data_nonlinear(q = q, p = p, n = n+n_test, m = 4, eff = eff, fixEff = fixEff)
    data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
    data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])

    fit <- SDForest(Y ~ ., data_train, cp = 0, mc.cores = mc.cores)
    fit2 <- ranger(Y ~ ., data_train, num.trees = 100, importance = 'impurity', mtry = floor(0.5 * ncol(data_train)))

    pred <- predict(fit, data_test)
    pred2 <- predict(fit2, data_test)$predictions

    mse <- (data$f_X[(n+1):(n+n_test)] - pred)
    mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)

    return(list(SDF = mse, ranger = mse2))
}