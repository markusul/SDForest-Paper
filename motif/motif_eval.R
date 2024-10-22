library(SDForest)
library(ggplot2)

load('motif/results/fits131.RData')

fit_dec
fit_plain

imp_dec <- fit_dec$var_importance
imp_dec <- imp_dec / max(imp_dec)

imp_plain <- fit_plain$var_importance
imp_plain <- imp_plain / max(imp_plain)

plot(imp_plain, imp_dec)

plot(sort(imp_dec), ylim = c(0, 1))
points(sort(imp_plain), col = 3)


path_dec <- regPath(fit_dec)
path_plain <- regPath(fit_plain)

stpath_dec <- stabilitySelection(fit_dec)
stpath_plain <- stabilitySelection(fit_plain)

plotOOB(path_dec)
plotOOB(path_plain)

plot(path_dec) + xlim(0, 0.03)
plot(path_plain) + xlim(0, 0.13)

plot(stpath_dec) + xlim(0, 0.03)
plot(stpath_plain) + xlim(0, 0.13)

fit_dec_pruned <- prune(fit_dec, path_dec$cp_min)
fit_plain_pruned <- prune(fit_plain, path_plain$cp_min)

imp_dec_pruned <- fit_dec_pruned$var_importance
imp_dec_pruned <- imp_dec_pruned / max(imp_dec_pruned)
names(imp_dec_pruned) <- 1:length(imp_dec_pruned)


imp_plain_pruned <- fit_plain_pruned$var_importance
imp_plain_pruned <- imp_plain_pruned / max(imp_plain_pruned)

plot(imp_plain_pruned, imp_dec_pruned)

plot(sort(imp_dec_pruned, decreasing = T))
plot(sort(imp_plain_pruned, decreasing = T))

which(imp_dec_pruned %in% sort(imp_dec_pruned, decreasing = T)[1:20])
imp_dec_pruned[41]

most_imp <- which.max(imp_dec_pruned)
most_imp

dep_dec <- partDependence(fit_dec_pruned, most_imp)
dep2_dec278 <- partDependence(fit_dec, 278, subSample = 500)
dep_dec278 <- partDependence(fit_dec_pruned, 278)
dep_plain278 <- partDependence(fit_plain_pruned, 278)
plot(dep_dec278)
plot(dep2_dec278)
plot(dep_plain278)

idx <- which(imp_dec_pruned %in% sort(imp_dec_pruned, decreasing = T)[1:20])
hist(apply(motifs, 2, var), breaks = 100)
idx




library(SDForest)
library(ggplot2)

load('motif/results/2fits131.RData')

path_dec <- regPath(fit_dec)
plotOOB(path_dec)

fit_dec <- prune(fit_dec, path_dec$cp_min)

path_plain <- regPath(fit_plain)
plotOOB(path_plain)

fit_plain <- prune(fit_plain, path_plain$cp_min)

imp_dec <- fit_plain$var_importance
names(imp_dec) <- 1:length(imp_dec)


plot(sort(fit_dec$var_importance, decreasing = T))
plot(sort(fit_plain$var_importance, decreasing = T))

sort(imp_dec, decreasing = T)[1:20]



#41        635         30        595         29        602        343        136         37 
#0.22495592 0.19590564 0.14524306 0.11424881 0.11000748 0.09047771 0.04226199 0.02669918 0.02541486 
#198        634        622        347        210        605        630        650        619 
#0.02539691 0.02132598 0.02001209 0.01725338 0.01547816 0.01434881 0.01385085 0.01339443 0.01157428 
#613         99 
#0.01150529 0.01102374


sort(imp_dec_pruned, decreasing = T)[1:20]
