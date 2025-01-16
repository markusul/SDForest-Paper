library(reticulate)
library(SDModels)
library(ranger)
np <- import("numpy")

# load data
npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]

# environments
table(interventions)

# set response variable as in DRIG
response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]

# load SDForest model trained on excluded interventions
load(file = 'cBench/results_rpe1_ex_sub/ENSG00000173812_sdf.Rdata')
fitsdf <- fromList(fitsdf)
fitsdf

# regularization path
path <- regPath(fitsdf)
plotOOB(path)
plot(path)

# stability selection
stabPath <- stabilitySelection(fitsdf)
plot(stabPath)

# estimate ranger model
fit_plain <- ranger(x = X[interventions == 'excluded', ], 
                    y = Y[interventions == 'excluded'], 
                    importance = 'impurity')
fit_plain

# variable importance
imp_plain <- fit_plain$variable.importance
imp_sdf <- fitsdf$var_importance

plot(imp_plain, imp_sdf)
plot(log(imp_plain), log(imp_sdf))

# predictions
pred_plain <- predict(fit_plain, X)
pred_plain <- pred_plain$predictions
pred_sdf <- predict(fitsdf, as.data.frame(X), mc.cores = 10)
save(pred_sdf, file = "cBench/results_rpe1_ex_sub/pred_sdf.RData")

plot(pred_plain, pred_sdf)
# add 90 degree line
abline(a = 0, b = 1, col = 'red')

residuals_plain <- pred_plain - Y
residuals_sdf <- pred_sdf - Y

plot(as.factor(interventions), residuals_plain)
plot(as.factor(interventions), residuals_sdf)

# MSE grouped per invtervention
env_MSE_plain <- tapply(residuals_plain^2, interventions, mean)
env_MSE_sdf <- tapply(residuals_sdf^2, interventions, mean)

plot(env_MSE_plain, env_MSE_sdf)
abline(a = 0, b = 1, col = 'red')
plot(log(env_MSE_plain), log(env_MSE_sdf))
abline(a = 0, b = 1, col = 'red')

plot(env_MSE_plain)
points(env_MSE_sdf, col = 'blue')

sort(env_MSE_plain - env_MSE_sdf)


sort(env_MSE_plain)

spec_envs <- which(names(env_MSE_plain) %in% c("excluded", "non-targeting", response))

par(mfrow = c(1, 2))
boxplot(env_MSE_plain[-spec_envs], ylim = c(0, 0.6))
grid()
boxplot(env_MSE_sdf[-spec_envs], ylim = c(0, 0.6))
grid()

env_MSE_plain[spec_envs]
env_MSE_sdf[spec_envs]

shift_order <- names(sort(env_MSE_plain[-spec_envs], decreasing = TRUE))

par(mfrow = c(1, 1))
plot(env_MSE_plain[shift_order])
points(env_MSE_sdf[shift_order], col = 'blue', pch = 3)

mean(env_MSE_plain[shift_order[1:50]])
mean(env_MSE_sdf[shift_order[1:50]])

plot((env_MSE_plain - env_MSE_sdf)[shift_order])
