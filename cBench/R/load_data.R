
library(SDForest)
library(reticulate)
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

fitsdf <- SDForest(x = X, y = Y, nTree = 100, gpu = TRUE, max_size = 1000)
plot(fitsdf$var_importance)
fitplain <- SDForest(x = X, y = Y, nTree = 100, gpu = TRUE, max_size = 1000, Q_type = "no_deconfounding")

plot(fitsdf$var_importance, fitplain$var_importance)

sort(fitsdf$var_importance, decreasing = T)[1:10]


path <- regPath(fitsdf)
