library(reticulate)
library(SDForest)
np <- import("numpy")

wass <- function(Fo, Y_int){
  
  Fi <- ecdf(Y_int)
  tempf <- function(x) abs(Fo(x) - Fi(x))
  val <- integrate(tempf,-Inf,Inf, subdivisions = 10000)$value
  return(val)
}

wass_xj <- function(Fo, Y, xj, interventions){
  idx <- interventions == xj
  wass(Fo, Y[idx])
}

wass2_xj <- function(Fo, Y, xj, interventions){
  idx <- interventions == xj
  wass2(Fo, Y[idx])
}

wass2 <- function(Y_obs, Y_int, p = 2){
  print(3)
  tempf <- function(x) abs(quantile(Y_obs, x, type = 5) - quantile(Y_int, x, type = 5))**p
  val <- integrate(tempf,0,1, subdivisions = 1000)$value **(1/p)
  return(val)
}

eff <- function(y_obs, y_int){
  (median(y_obs) - median(y_int))**2
}


eff_xj <- function(Y_obs, Y, xj, interventions){
  idx <- interventions == xj
  eff(Y_obs, 
      Y[idx])
}


npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]
sum(interventions == 'excluded')
sum(interventions == 'non-targeting')
response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
dim(X)

Y <- X[, response]
Y_obs <- Y[interventions == 'excluded']
Fo <- ecdf(Y_obs)

plot(svd(scale(X[interventions == 'excluded', ]))$d)
plot(svd(scale(X[interventions == 'excluded', ], scale = F))$d)
points(svd(scale(X[interventions == 'non-targeting', ], scale = F))$d, col = 2)

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

#load(file = 'cBench/results_rpe1/ENSG00000173812_sdf_ns.Rdata')
#fitsdf <- fromList(fitsdf)
#path <- regPath(fitsdf)
#plotOOB(path)
#path$cp_min
# 0.003

#fitsdf <- prune(fitsdf, 0.003)
#save(fitsdf, file = 'cBench/results_rpe1/ENSG00000173812_sdf_ns_pruned.Rdata')


load(file = 'cBench/results_rpe1_ex_sub/ENSG00000173812_sdf.Rdata')
load(file = 'cBench/results_rpe1_ex_sub/ENSG00000173812_plain.Rdata')

fitsdf <- fromList(fitsdf)
path <- regPath(fitsdf)



load(file = 'cBench/results_rpe1/ENSG00000173812_sdf_pruned.Rdata')
load(file = 'cBench/results_rpe1/ENSG00000173812_sdf_ns_pruned.Rdata')
load(file = 'cBench/results_rpe1/ENSG00000173812_plain_pruned_imp.Rdata')

imp_sdf <- fitsdf$var_importance
imp_plain <- fitsdf$var_importance

fit_plain <- ranger()


plot(imp_plain, imp_sdf)
plot(log(imp_plain), log(imp_sdf))

names(imp_sdf) <- fitsdf$var_names
names(imp_plain) <- fitsdf$var_names
imp_sdf_sort <- sort(imp_sdf, decreasing = T)
imp_plain_sort <- sort(imp_plain, decreasing = T)


top <- 20

sum(imp_sdf == 0)
plot(imp_plain, imp_sdf)
grid()
abline(h = imp_sdf_sort[top])
abline(v = imp_plain_sort[top])


plot(log(imp_plain), log(imp_sdf))
grid()
abline(h = log(imp_sdf_sort[top]))
abline(v = log(imp_plain_sort[top]))

var_names <- fitsdf$var_names
o_sdf <- var_names[imp_sdf >= imp_sdf_sort[top] & imp_plain < imp_plain_sort[top]]
o_plain <- var_names[imp_plain >= imp_plain_sort[top] & imp_sdf < imp_sdf_sort[top]]

wDist <- sapply(var_names, function(xj) wass_xj(Fo, Y, xj, interventions))
effDist <- sapply(var_names, function(xj) eff_xj(Y_obs, Y, xj, interventions))

plot(wDist, imp_plain)
plot(wDist, imp_sdf)

cor.test(wDist, imp_plain)
cor.test(wDist, imp_sdf)

plot(effDist, imp_plain)
plot(effDist, imp_sdf)

cor.test(effDist, imp_plain)
cor.test(effDist, imp_sdf)



mean(sapply(o_sdf, function(xj) wass_xj(Fo, Y, xj, interventions)))
mean(sapply(o_plain, function(xj) wass_xj(Fo, Y, xj, interventions)))

mean(sapply(o_sdf, function(xj) eff_xj(Y_obs, Y, xj, interventions)))
mean(sapply(o_plain, function(xj) eff_xj(Y_obs, Y, xj, interventions)))


Y[interventions == names(imp_plain_sort[3])]


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














transport::wasserstein(Y_obs, Y_int)

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]

sapply(names(sort(imp_sdf, decreasing = T)[1:20]), 
       function(xj) wass_xj(response, X, xj, interventions))

names(imp_plain) <- fitsdf$var_names
mean(sapply(names(sort(imp_plain, decreasing = T)[1:20]), 
           function(xj) cEff_xj(response, X, xj, interventions)), na.rm = T)
