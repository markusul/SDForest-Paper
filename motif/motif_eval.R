library(SDForest)
library(ggplot2)


load('motif/results/2fits131.RData')
fit_dec_ns <- copy(fit_dec)
load('motif/results/fits131.RData')

fit_dec
fit_dec_ns
fit_plain

path_dec <- regPath(fit_dec)
path_dec_ns <- regPath(fit_dec_ns)
path_plain <- regPath(fit_plain)

plotOOB(path_dec)
plotOOB(path_dec_ns)
plotOOB(path_plain)

plot(path_dec) + xlim(0, 0.03)
plot(path_dec_ns) + xlim(0, 0.03)
plot(path_plain) + xlim(0, 0.13)

fit_dec_pruned <- prune(fit_dec, path_dec$cp_min)
fit_dec_ns_pruned <- prune(fit_dec_ns, path_dec_ns$cp_min)
fit_plain_pruned <- prune(fit_plain, path_plain$cp_min)

imp_dec_pruned <- fit_dec_pruned$var_importance
names(imp_dec_pruned) <- 1:length(imp_dec_pruned)

imp_dec_ns_pruned <- fit_dec_ns_pruned$var_importance
names(imp_dec_ns_pruned) <- 1:length(imp_dec_ns_pruned)

imp_plain_pruned <- fit_plain_pruned$var_importance
names(imp_plain_pruned) <- 1:length(imp_plain_pruned)

plot(imp_plain_pruned / max(imp_plain_pruned), imp_dec_pruned /max(imp_dec_pruned), 
     ylab = 'deconfounded', xlab = 'plain')

plot(imp_dec_ns_pruned, imp_dec_pruned, 
     ylab = 'deconfounded', xlab = 'dec non scaled')

plot(imp_plain_pruned, imp_dec_ns_pruned, 
     ylab = 'deconfounded non scale', xlab = 'plain')

plot(sort(imp_dec_pruned, decreasing = T))
plot(sort(imp_plain_pruned, decreasing = T))

sort(imp_plain_pruned, decreasing = F)[1:40]


most_imp <- which.max(imp_dec_pruned)
most_imp

dep_dec <- partDependence(fit_dec_pruned, most_imp, subSample = 500)
dep_dec278 <- partDependence(fit_dec_pruned, 278, subSample = 500)
dep_dec_ns278 <- partDependence(fit_dec_ns_pruned, 278, subSample = 500)
dep_plain278 <- partDependence(fit_plain_pruned, 278, subSample = 500)
plot(dep_dec)
plot(dep_dec278)
plot(dep_dec_ns278)
plot(dep_plain278)



plot(fit_plain$oob_predictions, fit_plain$oob_predictions -fit_plain$Y)

plot(fit_dec$oob_predictions, fit_dec$oob_predictions -fit_dec$Y)


