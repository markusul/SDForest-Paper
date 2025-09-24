library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(reticulate)
library(SDModels)
library(ranger)
np <- import("numpy")

# singular values of cBench data
# load data
npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]

# set response variable as in DRIG
response <- "ENSG00000173812"

# only use observational data
env <- 'non-targeting'

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
Y <- X[interventions == env, response]
X <- X[interventions == env, -which(colnames(X) == response)]

X_dim <- dim(X)
lim <- (1+ sqrt(X_dim[2]/X_dim[1])) * sqrt(X_dim[1])

sig <- data.frame(d = svd(scale(X))$d)
sig$index <- 1:nrow(sig)

ggsig <- ggplot(sig, aes(y = d, x = index))+
  geom_point(size = 0.1) + 
  geom_abline(intercept = lim, slope = 0, lty = 2) + 
  ylab('singular values') + 
  xlim(1, max(sig$index) + 6)

ggsig
ggsave(filename = "simulation_study/figures/cBench_sig.jpeg", plot = ggsig, width = 3.5, height = 3)


# comparison of ranger to SDForest
load("cBench/semiSimResults/PlainSDF.RData")
load("cBench/semiSimResults/resRanger.RData")

par(mfrow = c(1, 3))
plot(imp_plain, imp_ranger)
plot(imp_sdf, imp_plain)
plot(imp_sdf, imp_ranger)


Plain <- imp_ranger
SDF <- imp_sdf
Plain <- Plain - min(Plain)
SDF <- SDF - min(SDF)
Plain <- Plain / max(Plain)
SDF <- SDF / max(SDF)

df_imp <- data.frame(Plain, SDF)
gg_imp <- ggplot(df_imp, aes(x = SDF, y = Plain)) +
  geom_point(size = 0.1) + 
  xlab('Variable importance SDForest') + 
  ylab('Variable importance ranger') + 
  geom_abline(lty = 2)
gg_imp
ggsave(filename = "simulation_study/figures/cBench_imp.jpeg", plot = gg_imp, width = 3.5, height = 3.5)

df_pred <- data.frame(preds_ranger, preds_sdf)
gg_pred <- ggplot(df_pred, aes(x = preds_sdf, y = preds_ranger)) + 
  geom_point(size = 0.1) + 
  xlab('Predictions of SDForest') + 
  ylab('Predictions of ranger') + 
  geom_abline(lty = 2) + xlim(2.4, 3) + ylim(2.4, 3)
gg_pred
ggsave(filename = "simulation_study/figures/cBench_pred.jpeg", plot = gg_pred, width = 3.5, height = 3.5)

library(gridExtra)
gg_imp <- gg_imp + annotate(geom = "text", label = 'a)', 
                            x = 0, 
                            y = 1)

gg_pred <- gg_pred + annotate(geom = "text", label = 'b)', 
                              x = 2.4, 
                              y = 3)

gg_comp <- grid.arrange(gg_imp, gg_pred, nrow = 1)
ggsave(filename = "simulation_study/figures/cBench_comp.jpeg", plot = gg_comp, width = 7, height = 4)
