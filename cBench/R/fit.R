library(reticulate)
library(SDForest)
np <- import("numpy")

npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
npz1$files
npz1$f[['var_names']]
interventions <- npz1$f[['interventions']]

response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
X <- X[interventions == 'non-targeting', ]
colnames(X) <- npz1$f[['var_names']]
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]

#RhpcBLASctl::omp_set_num_threads(1)

print("Fitting SDForest")
fitsdf <- SDForest(x = X, y = Y, nTree = 100, mc.cores = 100)
fitsdf <- toList(fitsdf)

print("Fitting SDForest with no deconfounding")
fitplain <- SDForest(x = X, y = Y, nTree = 100, mc.cores = 100,
                     Q_type = "no_deconfounding")
fitplain <- toList(fitplain)

print("Saving results")
save(fitplain, fitsdf, file = paste0('cBench/results_rpe1/', response, '.Rdata'))

