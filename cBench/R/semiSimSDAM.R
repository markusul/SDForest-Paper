set.seed(22)

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

# added linear dense confounding
p <- ncol(X)
n <- nrow(X)
q <- 1

Gamma <- matrix(rnorm(p*q), nrow = q)
delta <- matrix(rnorm(q), nrow = q)
H <- matrix(rnorm(q*n), nrow = n)

# without confounding
dat0 <- data.frame(Y, X)

# confounding strength
tau <- 1

Y_prime <- Y + H %*% delta * tau
X_prime <- X + H %*% Gamma * tau

dat <- data.frame(Y_prime, X_prime)

# fit SDAM
fit <- SDAM(Y ~ ., data = dat0, mc.cores = 20)
print(fit)
predsSDAM0 <- predict(fit, dat0)
fit <- NULL
fit <- SDAM(Y_prime ~ ., data = dat, mc.cores = 20)
print(fit)
predsSDAM <- predict(fit, dat)

print("save")
save(predsSDAM0, predsSDAM, file = paste0("cBench/semiSimResults/predSDAM.RData"))
