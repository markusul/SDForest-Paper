print('go')
library(anndata)
library(SDModels)
library(ranger)

npz1 <- read_h5ad("cBench/data/rpe1.h5ad")
print(npz1)
interventions <- npz1$obs[['gene']]

response <- "EIF1"

X <- npz1$X
X <- X[interventions == 'non-targeting', ]
colnames(X) <- npz1$var[['gene_name']]
print(colnames(X))
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]

#RhpcBLASctl::omp_set_num_threads(1)

print(dim(X))
print(apply(X, 2, function(x) {length(unique(x))}))

print("Fitting SDForest")
fitsdf <- SDForest(x = X, y = Y, nTree = 100, mc.cores = 100, max_size = 5000)
imp_sdf <- fitsdf$var_importance
preds_sdf <- predict(fitsdf, as.data.fram(npz1$X))

# estimate ranger model
fit_plain <- ranger(x = X, 
                    y = Y, 
                    importance = 'impurity')
preds_plain <- predict(fit_plain, npz1$X)
imp_plain <- fit_plain$variable.importance

Y <- npz1$X[, response]

print("Saving results")
save(imp_sdf, imp_plain, preds_plain, preds_sdf, interventions, Y, file = paste0('cBench/results_rpe1_full/', response, '.Rdata'))