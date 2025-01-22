library(reticulate)
library(SDModels)
library(ranger)
library(ggplot2)
library(ggsci)
np <- import("numpy")

# load data
npz1 <- np$load("cBench/data/dataset_rpe1_filtered.npz")
interventions <- npz1$f[['interventions']]

# environments
table(interventions)

# set response variable as in DRIG
response <- "ENSG00000173812"

X <- npz1$f[['expression_matrix']]
colnames(X) <- npz1$f[['var_names']]
Y <- X[, response]
X <- X[, -which(colnames(X) == response)]



# load SDForest model trained on excluded interventions
load(file = 'cBench/results_rpe1/ENSG00000173812_sdf.Rdata')
fitsdf <- fromList(fitsdf)
fitsdf

length(fitsdf$predictions)
sum(interventions == 'non-targeting')

# regularization path
#path <- regPath(fitsdf)
#plotOOB(path)
#plot(path)

# stability selection
#stabPath <- stabilitySelection(fitsdf)
#plot(stabPath)

# estimate ranger model
fit_plain <- ranger(x = X[interventions == 'non-targeting', ], 
                    y = Y[interventions == 'non-targeting'], 
                    importance = 'impurity')
fit_plain

# variable importance
imp_plain <- fit_plain$variable.importance
imp_sdf <- fitsdf$var_importance

plot(imp_plain, imp_sdf)
plot(log(imp_plain), log(imp_sdf))

# predictions
pred_plain <- predict(fit_plain, X)
pred_plain <- pred_plain$predictions
#pred_sdf <- predict(fitsdf, as.data.frame(X))
#save(pred_sdf, file = "cBench/results_rpe1/pred_sdf.RData")
load("cBench/results_rpe1/pred_sdf.RData")

plot(pred_plain[interventions == 'non-targeting'], pred_sdf[interventions == 'non-targeting'])
# add 90 degree line
abline(a = 0, b = 1, col = 'red')

plot(fit_plain$predictions, fitsdf$predictions)

residuals_plain <- pred_plain - Y
residuals_sdf <- pred_sdf - Y

env_id <- rep('blue', length(interventions))
env_id[interventions == 'excluded'] <- 'black'
env_id[interventions == 'non-targeting'] <- 'green'

env_name <- rep('env1', length(interventions))
env_name[interventions == 'excluded'] <- 'env2'
env_name[interventions == 'non-targeting'] <- 'env3'

par(mfrow = c(1, 2))
plot(pred_plain, residuals_plain, col = env_id)
points(pred_plain[interventions == 'non-targeting'], 
       residuals_plain[interventions == 'non-targeting'], col = 'green')
plot(pred_sdf, residuals_sdf, col = env_id)
points(pred_sdf[interventions == 'non-targeting'], 
       residuals_sdf[interventions == 'non-targeting'], col = 'green')


par(mfrow = c(1, 2))
plot(pred_plain, sqrt(abs(residuals_plain)), col = env_id)
points(pred_plain[interventions == 'non-targeting'], 
       sqrt(abs(residuals_plain[interventions == 'non-targeting'])), col = 'green')
plot(pred_sdf, sqrt(abs(residuals_sdf)), col = env_id)
points(pred_sdf[interventions == 'non-targeting'], 
       sqrt(abs(residuals_sdf[interventions == 'non-targeting'])), col = 'green')
title('scale-location')


model_analyse <- data.frame(prediction = c(pred_plain, pred_sdf), 
                            method = rep(c('ranger', 'SDF'), 
                                         each = length(pred_plain)), 
                            environment = as.factor(c(env_name, env_name)), 
                            residuals = c(residuals_plain, residuals_sdf), 
                            response = c(Y, Y))

# Calculate the size of each group
group_sizes <- table(model_analyse$environment)

model_analyse <- model_analyse[order(group_sizes[model_analyse$environment], decreasing = TRUE), ]

ggplot(model_analyse, aes(x = prediction, col = environment)) + 
  geom_point(aes(y = residuals)) + 
  facet_grid(~method) + 
  theme_bw() + 
  scale_color_tron()+ 
  ggtitle('Tukey-Anscombe')

ggplot(model_analyse, aes(x = prediction, col = environment)) + 
  geom_point(aes(y = sqrt(abs(residuals)))) + 
  facet_grid(~method) + 
  theme_bw() + 
  scale_color_tron() + 
  ggtitle('Scale-Location')



plot(as.factor(interventions), residuals_plain)
plot(as.factor(interventions), residuals_sdf)

# MSE grouped per invtervention
env_MSE_plain <- tapply(residuals_plain^2, interventions, mean)
env_MSE_sdf <- tapply(residuals_sdf^2, interventions, mean)
vars_Y <- tapply(Y, interventions, var)

plot(env_MSE_plain, env_MSE_sdf)
abline(a = 0, b = 1, col = 'red')
plot(log(env_MSE_plain), log(env_MSE_sdf))
abline(a = 0, b = 1, col = 'red')

plot(log(env_MSE_plain))
points(log(env_MSE_sdf), col = 'blue')

sort(env_MSE_plain - env_MSE_sdf)


sort(env_MSE_plain)

spec_envs <- which(names(env_MSE_plain) %in% c("excluded", "non-targeting", response))

par(mfrow = c(1, 2))
boxplot(env_MSE_plain[-spec_envs], ylim = c(0, 0.6))
grid()
boxplot(env_MSE_sdf[-spec_envs], ylim = c(0, 0.6))
grid()

env_MSE_plain[spec_envs]
env_MSE_sdf[spec_envs]

shift_order <- names(sort(env_MSE_plain[-spec_envs], decreasing = TRUE))

par(mfrow = c(1, 1))
plot(env_MSE_plain[shift_order])
points(env_MSE_sdf[shift_order], col = 'blue', pch = 3)

mean(env_MSE_plain[shift_order[1:50]])
mean(env_MSE_sdf[shift_order[1:50]])

plot((env_MSE_plain - env_MSE_sdf)[shift_order])
abline(h = 0, col = 'red')

boxplot(vars_Y[-spec_envs])


df_pred <- data.frame(pred_plain = pred_plain[interventions == 'non-targeting'], 
                      pred_sdf = pred_sdf[interventions == 'non-targeting'])

gg_pred <- ggplot(df_pred, aes(y = pred_sdf, x = pred_plain)) + 
  geom_point(size = 0.2) + theme_bw() + 
  ylab("prediction SDForest") + xlab("prediction ranger") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red")


df_imp <- data.frame(img_plain = imp_plain, 
                     imp_sdf = imp_sdf)

gg_imp <- ggplot(df_imp, aes(y = log(imp_sdf), x = log(imp_plain))) + 
  geom_point(size = 0.2) + theme_bw() + 
  ylab("log-importance SDForest") + xlab("log-importance ranger")

gg_imp
library(gridExtra)
gg_same <- grid.arrange(gg_imp, gg_pred, ncol = 2)
ggsave(filename = "simulation_study/figures/same.jpeg", plot = gg_same, width = 6, height = 3)

Y <- Y[interventions == 'non-targeting']
save(Y, file = "cBench/semiSimResults/Y.RData")
