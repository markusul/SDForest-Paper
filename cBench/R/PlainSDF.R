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

# classical random forest fit
fit_plain <- ranger(y = Y, x = X, importance = "impurity")
preds_plain <- fit_plain$predictions
imp_plain <- fit_plain$variable.importance

# deconfounded random forest fit
fit_sdf <- SDForest(y = Y, x = X, nTree = 500, 
                   max_size = 1000, mc.cores = 100)
preds_sdf <- fit_sdf$oob_predictions
imp_sdf <- fit_sdf$var_importance

print("save")
save(preds_plain, preds_sdf, imp_plain, imp_sdf, file = paste0("cBench/semiSimResults/PlainSDF.RData"))



