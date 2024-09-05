library(SDForest)
mc.cores <- 1

load("motif/data/MotifData.RData")

t <- 131
Y <- genes[, t]

fit_dec <- SDForest(x = motifs, y = Y, mc.cores = mc.cores)
fit_plain <- SDForest(x = motifs, y = Y, Q_type = 'no_deconfounding', mc.cores = mc.cores)

save(fit_dec, fit_plain, file = 'motif/results/fits.RData')
