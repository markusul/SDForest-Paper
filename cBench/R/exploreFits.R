library(SDForest)

#load(file = 'cBench/results_rpe1/ENSG00000173812_sdf.Rdata')
#load(file = 'cBench/results_rpe1/ENSG00000173812_plain.Rdata')

#fitplain <- fromList(fitplain)
#fitsdf <- fromList(fitsdf)

#path <- regPath(fitplain)
#path$cp_min
# 0.001

#fitplain <- prune(fitplain, 0.001)
#imp_plain <- fitplain$var_importance
#save(imp_plain, file = 'cBench/results_rpe1/ENSG00000173812_plain_pruned_imp.Rdata')

#path <- regPath(fitsdf)
#plotOOB(path)
#path$cp_min
# 0.002

#fitsdf <- prune(fitsdf, 0.002)
#save(fitsdf, file = 'cBench/results_rpe1/ENSG00000173812_sdf_pruned.Rdata')


load(file = 'cBench/results_rpe1/ENSG00000173812_sdf_pruned.Rdata')
load(file = 'cBench/results_rpe1/ENSG00000173812_plain_pruned_imp.Rdata')

imp_sdf <- fitsdf$var_importance

plot(imp_plain, imp_sdf)
length(imp_plain)

names(imp_sdf) <- fitsdf$var_names
sort(imp_sdf, decreasing = T)[1:5]

library(reticulate)
library(SDForest)
np <- import("numpy")

npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]

sort(imp_sdf, decreasing = T)[1:5]
response <- "ENSG00000173812"

xj <- names(sort(imp_sdf, decreasing = T))[3]
xj <- "ENSG00000108654"

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
X_int <- X[interventions == xj, ]
Y <- X_int[, response]
x <- X_int[, xj]

plot(Y ~ x)

X_obs <- X[interventions == 'non-targeting', ]
Y_obs <- X_obs[, response]
X_obs <- X_obs[, -which(colnames(X_obs) == response)]
x_obs <- X_obs[, xj]

X_rint <- X[interventions == response, ]
Y_rint <- X_rint[, response]
x_rint <- X_rint[, xj]

plot(Y_obs ~ x_obs)
points(Y ~ x, col = 'red')
points(Y_rint ~ x_rint, col = 'blue')

mean(Y_obs[x_obs == 0])
mean(Y[x == 0])

dep <- partDependence(fitsdf, xj, mc.cores = 5, subSample = 1000)
library(ggplot2)
plot(dep)

library(ranger)
fit <- ranger(x = X_obs, y = Y_obs, importance = 'impurity')
plot(imp_sdf, fit$variable.importance)
grid()
names(which(imp_sdf > 1e-4 & fit$variable.importance < 3.5))
names(which(imp_sdf > 9e-5 & fit$variable.importance < 2.5))
names(which(imp_sdf > 3e-4 & fit$variable.importance < 6))

pred <- predict(fit, X_int)$predictions

pred_sdf <- predict(fitsdf, data.frame(X_int))

plot(pred, pred_sdf)
plot(Y ~ x)
points(pred ~ x, col = 'red')
points(pred_sdf ~ x, col = 'blue')

mean(Y_obs[x_obs == 0])
mean(Y[x == 0])

mean(pred[x == 0])
mean(pred_sdf[x == 0])


plot(svd(scale(X_obs))$d)

cEff <- function(y_obs, y_int, xj_obs, xj_int){
  #(mean(y_obs[xj_obs == 0]) - mean(y_int[xj_int == 0]))**2
  (mean(y_obs) - mean(y_int))**2
}

cEff_xj <- function(response, X, xj, interventions){
  idx <- interventions == xj
  cEff(X[interventions == 'non-targeting', response], 
       X[idx, response], 
       X[interventions == 'non-targeting', xj], 
       X[idx, xj])
}

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]

mean(sapply(names(sort(imp_sdf, decreasing = T)[1:20]), 
       function(xj) cEff_xj(response, X, xj, interventions)), na.rm = T)

names(imp_plain) <- fitsdf$var_names
mean(sapply(names(sort(imp_plain, decreasing = T)[1:20]), 
           function(xj) cEff_xj(response, X, xj, interventions)), na.rm = T)
