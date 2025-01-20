library(ggplot2)

path <- "cBench/semiSimResults/"

load(paste0(path, "Y.RData"))

sims <- list.files(path)

res <- sapply(1:11, function(j){

  sim_j <- sims[grepl(paste0("predRob_", j, "_"), sims)]
  if(length(sim_j) < 11) return(NULL)

  Preds_plain <- matrix(NA, nrow = length(Y), ncol = length(sim_j))
  Preds_sdf <- matrix(NA, nrow = length(Y), ncol = length(sim_j))
  tau_seq <- numeric(length(sim_j))

  for (i in 1:length(sim_j)) {
    load(paste0(path, sim_j[i]))
    Preds_plain[, i] <- preds_plain
    Preds_sdf[, i] <- preds_sdf
    tau_seq[i] <- tau
  }

  rob_plain <- apply(Preds_plain, 2, 
                    function(pred) mean((pred - Preds_plain[, tau_seq == 0])**2))

  rob_sdf <- apply(Preds_sdf, 2, 
                    function(pred) mean((pred - Preds_sdf[, tau_seq == 0])**2))

  perf_plain <- apply(Preds_plain, 2, 
                    function(pred) mean((pred - Y)**2))

  perf_sdf <- apply(Preds_sdf, 2, 
                  function(pred) mean((pred - Y)**2))

  semiDat <- data.frame(perf = c(perf_plain, perf_sdf), 
                        rob = c(rob_plain, rob_sdf), 
                        tau = c(tau_seq, tau_seq), 
                        method = rep(c('ranger', 'SDF'), each = length(tau_seq)))
  semiDat
})

semiDat <- do.call(rbind, res)
semiDat$tau <- as.factor(semiDat$tau)

ggplot(semiDat, aes(x = tau, y = rob, col = method)) +
  geom_boxplot() + 
  ylim(0, 1) + 
  theme_bw() +
  ylab("model change")

ggplot(semiDat, aes(x = tau, y = perf, col = method)) +
  geom_boxplot() + 
  ylim(0, 1) +
  theme_bw() +
  ylab("performance change")
