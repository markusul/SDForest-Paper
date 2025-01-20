print('start')
library(reticulate)
library(SDModels)
np <- import("numpy")

npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]

response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
X <- X[interventions == 'excluded', ]
colnames(X) <- npz1$f[['var_names']]
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]

#RhpcBLASctl::omp_set_num_threads(1)

print("Fitting SDForest")
fitsdf <- SDForest(x = X, y = Y, nTree = 500, mc.cores = 100)
fitsdf <- toList(fitsdf)

print("Saving results")
save(fitsdf, file = paste0('cBench/results_rpe1_ex/', response, '_sdf.Rdata'))

#print("Fitting SDForest with no deconfounding")
#fitplain <- SDForest(x = X, y = Y, nTree = 500, mc.cores = 50,
#                     Q_type = "no_deconfounding")
#fitplain <- toList(fitplain)

#print("Saving results")
#save(fitplain, file = paste0('cBench/results_rpe1_ex/', response, '_plain.Rdata'))

