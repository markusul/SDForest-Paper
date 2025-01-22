library(ggplot2)
library(ggsci)

require(grid)
library(gridExtra)
grid_arrange_shared_legend <-
  function(...,
           ncol = length(list(...)),
           nrow = 1,
           position = c("bottom", "right"), 
           left = NULL, 
           bottom = NULL,
           vjust = 2.2) {
    
    plots <- list(...)
    position <- match.arg(position)
    g <-
      ggplotGrob(plots[[1]] + theme(legend.position = position))$grobs
    legend <- g[[which(sapply(g, function(x)
      x$name) == "guide-box")]]
    lheight <- sum(legend$height)
    lwidth <- sum(legend$width)
    gl <- lapply(plots, function(x)
      x + theme(legend.position = "none"))
    gl <- c(gl, ncol = ncol, nrow = nrow)
    
    combined <- switch(
      position,
      "bottom" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 1,
        heights = unit.c(unit(1, "npc") - lheight, lheight),
        left = grid.text(left, rot = 90, vjust = vjust), 
        bottom = grid.text(bottom)
      ),
      "right" = arrangeGrob(
        do.call(arrangeGrob, gl),
        legend,
        ncol = 2,
        widths = unit.c(unit(1, "npc") - lwidth, lwidth),
        left = grid.text(left, rot = 90, vjust = vjust),
        bottom = grid.text(bottom, vjust = -1.8)
      )
    )
    
    grid.newpage()
    grid.draw(combined)
    
    # return gtable invisibly
    invisible(combined)
    
  }

path <- "cBench/semiSimResults/"

load(paste0(path, "Y.RData"))

sims <- list.files(path)

res <- lapply(1:25, function(j){
  sim_j <- sims[grepl(paste0("predRob_", j, "_"), sims)]
  if(length(sim_j) < 11) {
    print(paste0(j, " not complete"))
    return(NULL)
  }

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

# diff ranger sdf at tau = 0
load(paste0(path, "predRob_1_1.RData"))
diff_rangerSDF <- mean((preds_plain - preds_sdf)**2)
diff_rangerSDF

gg_rob <- ggplot(semiDat, aes(x = tau, y = rob, col = method)) +
  geom_boxplot(outlier.size = 0.4) + 
  theme_bw() +
  ylab("model change") + 
  scale_color_tron() +
  geom_boxplot(aes(x = "0", y = diff_rangerSDF, col = 'SDF - ranger')) +
  theme(legend.title = element_blank())
gg_rob
ggsave(filename = "simulation_study/figures/SemiSim_rob.jpeg", 
       plot = gg_rob, width = 6, height = 4)

gg_rob_log <- ggplot(semiDat, aes(x = tau, y = log(rob), col = method)) +
  geom_boxplot(outlier.size = 0.4) + 
  theme_bw() +
  ylab("log model change") + 
  scale_color_tron() +
  geom_boxplot(aes(x = "0", y = log(diff_rangerSDF), col = 'SDF - ranger')) +
  theme(legend.title = element_blank())
gg_rob_log
ggsave(filename = "simulation_study/figures/SemiSim_rob_log.jpeg", 
       plot = gg_rob_log, width = 6, height = 4)

gg_perf <- ggplot(semiDat, aes(x = tau, y = perf, col = method)) +
  geom_boxplot(outlier.size = 0.4) + 
  ylim(0, 1) + 
  theme_bw() +
  ylab("performance") + 
  scale_color_tron() +
  theme(legend.title = element_blank())

gg_semiSim <- grid_arrange_shared_legend(gg_rob, gg_rob_log, ncol = 2)
ggsave(filename = "simulation_study/figures/SemiSim.jpeg", 
       plot = gg_semiSim, width = 10, height = 4)



# performance with tau = 0
mean(semiDat[semiDat$tau == 0 & semiDat$method == "ranger", "perf"])
mean(semiDat[semiDat$tau == 0 & semiDat$method == "SDF", "perf"])
var(Y)

1 - mean(semiDat[semiDat$tau == 0 & semiDat$method == "ranger", "perf"]) / var(Y)
1 - mean(semiDat[semiDat$tau == 0 & semiDat$method == "SDF", "perf"]) / var(Y)


load(paste0(path, "predRob_1_1.RData"))

plot(preds_plain)

