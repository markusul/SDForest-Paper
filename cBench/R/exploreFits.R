library(SDForest)

load(file = 'cBench/results_rpe1/ENSG00000173812_sdf.Rdata')
fitplain <- fromList(fitplain)
fitsdf <- fromList(fitsdf)

plot(fitsdf$var_importance, fitplain$var_importance)


path <- regPath(fitsdf)
