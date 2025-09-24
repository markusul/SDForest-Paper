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

# only use observational data
env <- 'non-targeting'

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
Y <- X[interventions == env, response]
X <- X[interventions == env, -which(colnames(X) == response)]

start.time <- Sys.time()

# classical random forest fit
fit_ranger <- ranger(y = Y, x = X, importance = "impurity", 
                     mtry = floor(0.5 * ncol(X)),
                     sample.fraction = 1000/nrow(X))
preds_ranger <- fit_ranger$predictions
imp_ranger <- fit_ranger$variable.importance
save(preds_ranger, imp_ranger, file = paste0("cBench/semiSimResults/resRanger.RData"))
print(Sys.time() - start.time)

fit <- SDForest(y = Y, x = X, nTree = 500, 
                   max_size = 1000, mc.cores = 100, 
                   Q_type = "no_deconfounding")
preds_plain <- fit$oob_predictions
imp_plain <- fit$var_importance

fit <- NULL
gc()
print(Sys.time() - start.time)

# deconfounded random forest fit
fit <- SDForest(y = Y, x = X, nTree = 500, 
                   max_size = 1000, mc.cores = 100)
preds_sdf <- fit$oob_predictions
imp_sdf <- fit$var_importance

print(Sys.time() - start.time)
print("save")
save(preds_plain, preds_sdf, imp_plain, imp_sdf, file = paste0("cBench/semiSimResults/PlainSDF.RData"))




