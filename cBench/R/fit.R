library(reticulate)
library(SDForest)
np <- import("numpy")

npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
npz1$files
npz1$f[['var_names']]
interventions <- npz1$f[['interventions']]

dim(npz1$f[['expression_matrix']])


response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
X <- X[interventions == 'non-targeting', ]
colnames(X) <- npz1$f[['var_names']]
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]

#RhpcBLASctl::omp_set_num_threads(1)

fitsdf <- SDForest(x = X, y = Y, nTree = 100, mc.cores = 100)
fitsdf <- toList(fitsdf)

fitplain <- SDForest(x = X, y = Y, nTree = 100, mc.cores = 100,
                     Q_type = "no_deconfounding")
fitplain <- toList(fitplain)

save(fitplain, fitsdf, file = 'cBench/results/fits.Rdata')

