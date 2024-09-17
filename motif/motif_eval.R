library(SDForest)

load('motif/results/fits.RData')

fit_dec
fit_plain

imp_dec <- fit_dec$var_importance
imp_dec <- imp_dec / max(imp_dec)

imp_plain <- fit_plain$var_importance
imp_plain <- imp_plain / max(imp_plain)

plot(imp_plain, imp_dec)

plot(sort(imp_dec))
points(sort(imp_plain), col = 3)


path_dec <- regPath(fit_dec)
path_plain <- regPath(fit_plain)

stpath_dec <- stabilitySelection(fit_dec)
stpath_plain <- stabilitySelection(fit_plain)


plot(path_dec)
plot(path_plain)

plot(stpath_dec)
plot(stpath_plain)

