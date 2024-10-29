library(SDForest)

load(file = 'cBench/results/fits.Rdata')
fitplain <- fromList(fitplain)
fitsdf <- fromList(fitsdf)

plot(fitsdf$var_importance, fitplain$var_importance)
