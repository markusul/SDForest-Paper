args = commandArgs(trailingOnly = TRUE)
set.seed(as.numeric(args[2]))

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

env <- 'non-targeting'

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
Y <- X[interventions == env, response]
X <- X[interventions == env, -which(colnames(X) == response)]

p <- ncol(X)
n <- nrow(X)
q <- 1

Gamma <- matrix(rnorm(p*q), nrow = q)
delta <- matrix(rnorm(q), nrow = q)
H <- matrix(rnorm(q*n), nrow = n)

tau_seq <- seq(0, 1, 0.1)
tau <- tau_seq[as.numeric(args[1])]

Y_prime <- Y + H %*% delta * tau
X_prime <- X + H %*% Gamma * tau
fitTau <- ranger(y = Y_prime, x = X_prime)
preds_plain <- fitTau$predictions

fitTau <- SDForest(y = Y_prime, x = X_prime, nTree = 500, 
                   max_size = 1000, mc.cores = 100)
preds_sdf <- fitTau$oob_predictions

#preds_plain <- sapply(tau_seq, function(tau){
#  Y_prime <- Y + H %*% delta * tau
#  X_prime <- X + H %*% Gamma * tau
#  fitTau <- ranger(y = Y_prime, x = X_prime)
#  fitTau$predictions
#})

#preds_sdf <- sapply(tau_seq, function(tau){
#  Y_prime <- Y + H %*% delta * tau
#  X_prime <- X + H %*% Gamma * tau
#  fitTau <- SDForest(y = Y_prime, x = X_prime, nTree = 500, 
#                     max_size = 1000, mc.cores = 100)
#  fitTau$oob_predictions
#})

print("save")
save(preds_plain, preds_sdf, tau, file = paste0("cBench/semiSimResults/predRob_", args[2], "_", args[1], ".RData"))


