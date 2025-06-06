library(SDModels)
mc.cores <- 100

p <- 500
n <- 1000
q <- 20

n_test <- 500

set.seed(2024)
data <- simulate_data_nonlinear(q, p, n + n_test, 4)
data_test <- data
data_test$Y <- data_test$Y[(n+1):(n+n_test)]
data_test$X <- data_test$X[(n+1):(n+n_test),]
data_test$f_X <- data_test$f_X[(n+1):(n+n_test)]

data$X <- data$X[1:n,]
data$Y <- matrix(data$Y[1:n])
data$f_X <- data$f_X[1:n]

start <- Sys.time()
fit <- SDForest(x = data$X, y = data$Y, cp = 0, mc.cores = mc.cores)
end <- Sys.time()
print(end - start)

print('fit done')

reg_path <- regPath(fit)
print('reg_path done')

stable_path <- stabilitySelection(fit)
print('stable_path done')

dep_1 <- partDependence(fit, data$j[1], mc.cores = mc.cores)
dep_2 <- partDependence(fit, data$j[2], mc.cores = mc.cores)
dep_3 <- partDependence(fit, data$j[3], mc.cores = mc.cores)
dep_4 <- partDependence(fit, data$j[4], mc.cores = mc.cores)
print('dep done')

save(data, data_test, fit, reg_path, stable_path, dep_1, dep_2, dep_3, dep_4, 
    file = "simulation_study/results/default_szenario.RData")
print('save done')
