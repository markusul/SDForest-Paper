library(SDForest)
mc.cores <- 5

load("motif/data/MotifData.RData")

t <- 131
Y <- genes[, t]

fit_plain <- SDForest(x = motifs, y = Y, Q_type = 'no_deconfounding', mc.cores = mc.cores, Q_scale = F)
print(fit_plain)
fit_dec <- SDForest(x = motifs, y = Y, mc.cores = mc.cores, Q_scale = F)
print(fit_dec)

save(fit_dec, fit_plain, file = paste0('motif/results/2fits', as.character(t) ,'.RData'))
