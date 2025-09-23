require(grid)
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

true_function <- function(beta, js){
    res <- list(beta = beta, js = js)
    class(res) <- 'true_function'
    return(res)
}

predict.true_function <- function(object, newdata){
    f_X <- apply(newdata, 1, function(x) f_four(x, object$beta, object$js))
    return(f_X)
}

plotDep <- function(object, n_examples = 19){
  ggdep <- ggplot2::ggplot()
  preds <- object$preds
  x_seq <- object$x_seq
  
  sample_examples <- sample(1:dim(preds)[2], n_examples)
  for(i in sample_examples){
      pred_data <- data.frame(x = x_seq, y = preds[, i])
      ggdep <- ggdep + ggplot2::geom_line(data = pred_data, ggplot2::aes(x = x, y = y), col = 'grey')
  }

  ggdep <- ggdep + ggplot2::geom_line(data = data.frame(x = x_seq, y = object$preds_mean), 
                    ggplot2::aes(x = x, y = y), col = '#08cbba', linewidth = 1.5)
  ggdep <- ggdep + ggplot2::geom_point(data = data.frame(x = object$xj, y = -5), 
                    ggplot2::aes(x = x, y = y), col = 'black', size = 1,shape = 108)
  ggdep <- ggdep + ggplot2::ylab('')
  ggdep <- ggdep + ggplot2::xlim(quantile(object$xj, 0.05), quantile(object$xj, 0.95))
  if(is.character(object$j)){
    ggdep <- ggdep + ggplot2::xlab(object$j)
  }else{
    ggdep <- ggdep + ggplot2::xlab(paste('x', object$j, sep = ''))
  }
  ggdep + ggplot2::ylim(-2.2, 2)
}

ranger_fun <- function(object){
    res <- list(model = object)
    class(res) <- 'ranger_fun'
    return(res)
}
predict.ranger_fun <- function(object, newdata){
    return(predict(object$model, newdata)$predictions)
}

library(SDModels)
library(gridExtra)
library(ggplot2)
theme_set(theme_bw(base_size = 14))
library(ranger)
library(ggsci)
library(ggpubr)
library(tidyr)

##### motivation of loss #####

set.seed(22)
dat <- simulate_data_nonlinear(q = 20, p = 500, n = 1000, m = 4) 
Y <- dat$Y
f <- dat$f_X
Q <- get_Q(dat$X, 'trim')
QY <- Q %*% Y
Qf <- Q %*% f

df <- data.frame(f, Y, QY, Qf)

plain <- ggplot(df, aes(y = f, x = Y)) + 
  geom_point(size = 0.2) + 
  ylab("f(X)") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  ylim(-5, 5) + xlim(-20, 20)

transformed <- ggplot(df, aes(y = Qf, x = QY)) + 
  geom_point(size = 0.2) + 
  ylab("Qf(X)") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  ylim(-5, 5) + xlim(-20, 20)

pt <- grid.arrange(plain, transformed, ncol = 2)
ggsave(filename = "simulation_study/figures/pt.jpeg", plot = pt, width = 7, height = 3)

# show shrinkage of singular values
set.seed(22)
q = 20
dat <- simulate_data_nonlinear(q = q, p = 100, n = 1000, m = 4) 
X <- scale(dat$X)

d <- svd(X)$d

QTrim <- get_Q(X, type = 'trim', scaling = F)
dTrim <- svd(QTrim %*% X)$d

QPCA <- get_Q(X, type = 'pca', scaling = F, q_hat = q)
dPCA <- svd(QPCA %*% X)$d
dPCA <- c(rep(0, q), dPCA[1:(length(dPCA)-q)])

sigma_dat <- data.frame(sigma = c(d, dTrim, dPCA), 
                        method = rep(c('no', 'trim', 'pca'), each = length(d)), 
                        index = rep(1:length(d), 3))

ggsig <- ggplot(sigma_dat, aes(y = sigma, x = index, col = method, shape = method)) + 
  geom_point(size = 0.7) + 
  ylab('transformed singular values')
  
ggsig
ggsave(filename = "simulation_study/figures/sig.jpeg", plot = ggsig, width = 7, height = 3)


##### default experiment #####

load('simulation_study/results/default_szenario.RData')
plotOOB(reg_path)

set.seed(2024)

fit2 <- ranger(x = data.frame(data$X), y = data$Y, num.trees = 100, 
  importance = 'impurity', mtry = floor(0.5 * ncol(data$X)))

true_f <- true_function(data$beta, data$j)

dep_f_1 <- partDependence(true_f, data$j[1], data$X)
dep_f_2 <- partDependence(true_f, data$j[2], data$X)
dep_f_3 <- partDependence(true_f, data$j[3], data$X)
dep_f_4 <- partDependence(true_f, data$j[4], data$X)

imp_1 <- fit$var_importance / max(fit$var_importance)
imp_2 <- fit2$variable.importance / max(fit2$variable.importance)
true_imp <- rep('spurious     ', length(imp_1))
true_imp[data$j] <- 'causal'

#order of causal parents in the model
fileConn<-file("simulation_study/figures/causal_order.txt")

# SDF
sdf_order <- paste(c('SDF:', which(true_imp[order(imp_1, decreasing = T)] == 'causal')), 
                   collapse = ' ')
#ranger
ranger_order <- paste(c('ranger:', which(true_imp[order(imp_2, decreasing = T)] == 'causal')), 
                   collapse = ' ')
writeLines(c(sdf_order, ranger_order), fileConn)
close(fileConn)

imp_data <- data.frame(SDF = imp_1, ranger = imp_2, Covariates = as.factor(true_imp))



ggimp <- ggplot(imp_data, aes(x = SDF, y = ranger, 
                              col = Covariates, shape = Covariates)) + 
  geom_point(size = 0.7) + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Normalized to [0, 1]') + 
  theme(legend.title = element_blank())

ggimp_log <- ggplot(imp_data, aes(x = log(SDF), y = log(ranger), 
                                  col = Covariates, shape = Covariates)) + 
  geom_point(size = 0.7) + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Logarithmic scale') + 
  theme(legend.title = element_blank())

gg_imp <- grid_arrange_shared_legend(ggimp, ggimp_log, position = 'right',
                           left = 'Variable importance ranger', 
                           bottom = 'Variable importance SDForest')

ggsave(filename = "simulation_study/figures/imp.jpeg", plot = gg_imp, width = 8, height = 4)

ggimp_log <- ggplot(imp_data, aes(x = log(SDF), y = log(ranger), 
                                  col = Covariates, shape = Covariates)) + 
  geom_point(size = 0.7) + xlab('log Variable importance SDForest') + 
  ylab('log Variable importance ranger') + scale_color_tron() + 
  theme(legend.title = element_blank(), legend.position = "inside",
        legend.position.inside=c(0.8,0.8)) + 
  scale_shape_manual(values=c(3, 16))
ggimp_log
ggsave(filename = "simulation_study/figures/imp_log.jpeg", plot = ggimp_log, width = 5, height = 4)

ggdep1 <- plotDep(dep_1) + 
  geom_line(aes(x = dep_f_1$x_seq, y = dep_f_1$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep2 <- plotDep(dep_2) + 
  geom_line(aes(x = dep_f_2$x_seq, y = dep_f_2$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep3 <- plotDep(dep_3) +
  geom_line(aes(x = dep_f_3$x_seq, y = dep_f_3$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep4 <- plotDep(dep_4) +
  geom_line(aes(x = dep_f_4$x_seq, y = dep_f_4$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

gg_cond_rf <- grid_arrange_shared_legend(ggdep1, ggdep2, ggdep3, ggdep4,
                           ncol = 2, nrow = 2, left = 'f(X)')
ggsave(filename = "simulation_study/figures/cond_rf.jpeg", 
       plot = gg_cond_rf, width = 7, height = 6)

ranger_fit <- ranger_fun(fit2)

dep_r_1 <- partDependence(ranger_fit, data$j[1], data.frame(data$X))
dep_r_2 <- partDependence(ranger_fit, data$j[2], data.frame(data$X))
dep_r_3 <- partDependence(ranger_fit, data$j[3], data.frame(data$X))
dep_r_4 <- partDependence(ranger_fit, data$j[4], data.frame(data$X))

ggdep1_r <- plotDep(dep_r_1) + 
  geom_line(aes(x = dep_f_1$x_seq, y = dep_f_1$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep2_r <- plotDep(dep_r_2) +
  geom_line(aes(x = dep_f_2$x_seq, y = dep_f_2$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep3_r <- plotDep(dep_r_3) + 
  geom_line(aes(x = dep_f_3$x_seq, y = dep_f_3$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

ggdep4_r <- plotDep(dep_r_4) +
  geom_line(aes(x = dep_f_4$x_seq, y = dep_f_4$preds_mean, col = 'red'), 
            linetype = 2, linewidth = 0.2) + 
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c("red"), labels = c("True Function"))

gg_cond_r <- grid_arrange_shared_legend(ggdep1_r, ggdep2_r, ggdep3_r, ggdep4_r,
                                         ncol = 2, nrow = 2, left = 'f(X)')
ggsave(filename = "simulation_study/figures/cond_r.jpeg", 
       plot = gg_cond_r, width = 3.5, height = 4)

gg_regpath <- ggplot()
for(i in 1:ncol(reg_path$varImp_path)){
  gg_regpath <- gg_regpath + geom_line(data = data.frame(x = reg_path$cp, 
    y = reg_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey', 
    linewidth = if(i %in% data$j) 0.7 else 0.2)
}
gg_regpath <- gg_regpath + xlab('') + 
  ylab('Variable importance') + ggtitle('Variable importance path') +
  xlim(0, 0.4) + geom_line(aes(x = 1, y = 0, col = '#d11010'), linewidth = 0.7) +
  geom_line(aes(x = 1, y = 0, col = 'grey'), linewidth = 0.2) +
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c('#d11010', 'grey'), 
                              labels = c("causal \nparents     ", 'other \nvariables'))

gg_stablepath <- ggplot()
for(i in 1:ncol(stable_path$varImp_path)){
  gg_stablepath <- gg_stablepath + geom_line(data = data.frame(x = stable_path$cp, 
    y = stable_path$varImp_path[, i]), aes(x = x, y = y), 
    col = if(i %in% data$j)'#d11010' else 'grey', 
    linewidth = if(i %in% data$j) 0.7 else 0.2)
}
gg_stablepath <- gg_stablepath + xlab('') + 
  ylab(expression(Pi)) + ggtitle('Stability selection path') +
  xlim(0, 0.4) + geom_line(aes(x = 1, y = 0, col = '#d11010'), linewidth = 0.7) +
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c('#d11010', 'grey'), 
                              labels = c("causal \nparents    ", 'other \nvariables'))

gg_paths <- grid_arrange_shared_legend(gg_regpath, gg_stablepath, ncol = 2, 
                                       bottom = 'Regularization: cp', 
                                       position = 'right')

ggsave(filename = "simulation_study/figures/paths.jpeg", 
       plot = gg_paths, width = 9, height = 5)

##### Performance depending on the dimensions #####
error_name <- expression("||"*f^0*(x[test]) - hat(f)(x[test])*"||"[2]^2 / n[test])
annot_x_shift <- 0.23
annot_y_shift <- 1

agg_fun <- function(x){
  mean(x**2)
}

load_perf <- function(path, agg_fun){
  load(path)
  perf <- lapply(perf, function(n) cbind(t(sapply(n, function(x)sapply(x, agg_fun))), seq = seq_))

  perf <- do.call(rbind, perf)
  perf <- data.frame(perf)
  perf$seq <- as.factor(perf$seq)

  perf <- gather(perf, 'method', error, -seq)
  perf
}

#### n
files <- list.files('simulation_study/results/perf_n')
length(files)
perf_n <- lapply(paste0('simulation_study/results/perf_n/', files), 
  load_perf, agg_fun = agg_fun)

perf_n <- do.call(rbind, perf_n)

gg_n <- ggplot(perf_n, aes(x = seq, y = error, color = method)) + 
  geom_boxplot(outlier.size = 0.4) + xlab('Number of training samples') + 
  ylab('') + scale_color_tron() + theme(legend.title=element_blank())
  
gg_n <- gg_n + annotate(geom = "text", label = 'a)', 
                x = ggplot_build(gg_n)$layout$panel_scales_x[[1]]$range_c$range[[1]] + annot_x_shift, 
                y = ggplot_build(gg_n)$layout$panel_scales_y[[1]]$range$range[[2]] - annot_y_shift)
gg_n

#### p
files <- list.files('simulation_study/results/perf_p')
length(files)

perf_p <- lapply(paste0('simulation_study/results/perf_p/', files), 
  load_perf, agg_fun = agg_fun)

perf_p <- do.call(rbind, perf_p)

gg_p <- ggplot(perf_p, aes(x = seq, y = error, color = method)) + 
  geom_boxplot(outlier.size = 0.4) + xlab('Number of covariates') + 
  ylab('') + scale_color_tron() + theme(legend.title=element_blank())

gg_p <- gg_p + annotate(geom = "text", label = 'b)', 
                        x = ggplot_build(gg_p)$layout$panel_scales_x[[1]]$range_c$range[[1]] + 0.5 + annot_x_shift, 
                        y = ggplot_build(gg_p)$layout$panel_scales_y[[1]]$range$range[[2]] - annot_y_shift)
gg_p

#### q
files <- list.files('simulation_study/results/perf_q')
length(files)

perf_q <- lapply(paste0('simulation_study/results/perf_q/', files), 
  load_perf, agg_fun = agg_fun)

perf_q <- do.call(rbind, perf_q)

gg_q <- ggplot(perf_q, aes(x = seq, y = error, color = method)) + 
  geom_boxplot(outlier.size = 0.4) + xlab('Number of confounders') + 
  ylab('') + scale_color_tron() + theme(legend.title=element_blank())

gg_q <- gg_q + annotate(geom = "text", label = 'c)', 
                        x = ggplot_build(gg_q)$layout$panel_scales_x[[1]]$range_c$range[[1]] + 0.05 + annot_x_shift, 
                        y = ggplot_build(gg_q)$layout$panel_scales_y[[1]]$range$range[[2]] - annot_y_shift)
gg_q

#### Density assumption ####
files <- list.files('simulation_study/results/perf_eff')
length(files)
files <- paste0(1:10, ".RData")

perf_eff <- lapply(paste0('simulation_study/results/perf_eff/', files), 
  load_perf, agg_fun = agg_fun)

perf_eff <- do.call(rbind, perf_eff)

gg_eff <- ggplot(perf_eff, aes(x = seq, y = error, color = method)) + 
  geom_boxplot(outlier.size = 0.4) + xlab("Number of affected covariates") + 
  ylab('') + scale_color_tron() + theme(legend.title=element_blank())

gg_eff <- gg_eff + annotate(geom = "text", label = 'd)', 
                        x = ggplot_build(gg_eff)$layout$panel_scales_x[[1]]$range_c$range[[1]] + 0.5 + annot_x_shift, 
                        y = ggplot_build(gg_eff)$layout$panel_scales_y[[1]]$range$range[[2]] - annot_y_shift)
gg_eff

gg_dims <- grid_arrange_shared_legend(gg_n, gg_p, gg_q, gg_eff, 
                           ncol = 2, nrow = 2, 
                           left = error_name, vjust = 1.5)
gg_dims

ggsave(filename = "simulation_study/figures/dims1.jpeg", 
       plot = gg_dims, width = 8, height = 6)

dim_names <- c('Number of training samples', 'Number of covariates', 
               'Number of confounders', "Number of affected covariates")

perf_n$dim <- dim_names[1]
perf_p$dim <- dim_names[2]
perf_q$dim <- dim_names[3]
perf_eff$dim <- dim_names[4]

perf_dim <- rbind(perf_n, perf_p, perf_q, perf_eff)
perf_dim$seq <- factor(perf_dim$seq, order = TRUE, 
                       levels = as.character(sort(as.numeric(levels(perf_dim$seq)))))
perf_dim$dim <- factor(perf_dim$dim, ordered = T, levels = dim_names)

gg_dims2 <- ggplot(perf_dim, aes(x = seq, y = error, color = method)) + 
  geom_boxplot(outlier.size = 0.4) + xlab("") + 
  ylab(error_name) + scale_color_tron() + facet_wrap(~dim, ncol = 2, scales="free") + 
  theme(legend.position = 'bottom', legend.title = element_blank())

dat_text <- data.frame(
  label = factor(c("a)", "b)", "c)", "d)"), ordered = T),
  dim   = factor(dim_names, ordered = T),
  method = 'ranger'
)

gg_dims2 <- gg_dims2 + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = Inf, 
                label = label, 
                vjust = 1.5, hjust = -0.6), color = 'black'
)
gg_dims2

ggsave(filename = "simulation_study/figures/dims2.jpeg", 
       plot = gg_dims2, width = 8, height = 6)

# nonlinear confounding with fourier
load('simulation_study/results/nonlin_confounding_2.RData')

gg_sing <- ggplot(sing, aes(x = i, y = d)) + 
  geom_point(size = 0.7) +
  ylab("singular values") + 
  xlab("index")
gg_sing
ggsave(filename = "simulation_study/figures_nl/sing2.jpeg", plot = gg_sing, width = 7, height = 3)


plain <- ggplot(df, aes(y = f, x = Y)) + 
  geom_point(size = 0.2) +
  ylab("f(X)") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  ylim(-1.5, 0.8) + xlim(-7, 7)

transformed <- ggplot(df, aes(y = Qf, x = QY)) + 
  geom_point(size = 0.2) +
  ylab("Qf(X)") + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red") + 
  ylim(-1.5, 0.8) + xlim(-7, 7)

pt <- grid.arrange(plain, transformed, ncol = 2)
ggsave(filename = "simulation_study/figures_nl/pt_nl2.jpeg", plot = pt, width = 7, height = 3)


# Comparison of Variable importance
imp_1 <- fit$var_importance / max(fit$var_importance)
imp_2 <- fit2$var_importance / max(fit2$var_importance)
true_imp <- rep('spurious     ', length(imp_1))
true_imp[js] <- 'causal'

imp_data <- data.frame(SDF = imp_1, ranger = imp_2, Covariates = as.factor(true_imp))

ggimp <- ggplot(imp_data, aes(x = SDF, y = ranger, col = Covariates)) + 
  geom_point(size = 0.5) + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Normalized to [0, 1]') + 
  theme(legend.title = element_blank())

ggimp_log <- ggplot(imp_data, aes(x = log(SDF), y = log(ranger), col = Covariates)) + 
  geom_point(size = 0.5) + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Logarithmic scale') + 
  theme(legend.title = element_blank())

gg_imp <- grid_arrange_shared_legend(ggimp, ggimp_log, position = 'right',
                                     left = 'Variable importance ranger', 
                                     bottom = 'Variable importance SDForest')
ggsave(filename = "simulation_study/figures_nl/imp_nl2.jpeg", plot = gg_imp, width = 8, height = 4)

set.seed(22)
gg_dep <- plot(dep) + geom_point(aes(x = X[, js[1]], y = Y), size = 0.2)

sample_examples <- sample(1:ncol(dep2$preds), 19)
for(i in sample_examples){
  pred_data <- data.frame(x = dep2$x_seq, y = dep2$preds[, i])
  gg_dep <- gg_dep + ggplot2::geom_line(data = pred_data, 
                                        ggplot2::aes(x = x, y = y), col = 'lightblue')
}
gg_dep <- gg_dep + 
  geom_line(aes(x = dep2$x_seq, y = dep2$preds_mean, col = 'rf')) + 
  geom_point(aes(x = X[, js[1]], y = f_X, col = 'true'), size = 0.2) +   
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c(true = "red", rf = "blue"), 
                              labels = c(true = "True Function", 
                                         rf = "no deconfounding")) +
  theme(legend.position = 'bottom')
gg_dep

set.seed(22)
gg_dep <- plot(dep) + geom_point(aes(x = X[, js[1]], y = Y), size = 0.2)
gg_dep <- gg_dep + 
  geom_line(aes(x = X[, js[1]], y = f_X, col = 'true'), 
            linetype = 2, linewidth = 0.2) +   
  ggplot2::labs(col = "") + 
  ggplot2::scale_color_manual(values = c(true = "red"), 
                              labels = c(true = "True Function")) +
  theme(legend.position = 'bottom')
gg_dep
ggsave(filename = "simulation_study/figures_nl/dep_nl2.jpeg", plot = gg_dep, width = 7, height = 4)


gg_perf <- ggplot(perf_g, aes(y = performance, x = method, color = method)) +
  geom_boxplot(outlier.size = 0.4) +
  scale_fill_tron() + ylab(error_name) + xlab('') +
  theme(legend.position = 'None')

gg_perf
ggsave(filename = "simulation_study/figures_nl/perf_nl2.jpeg", plot = gg_perf, width = 3.5, height = 3)


#### fast comparison ####
load('simulation_study/results/fast.RData')

gg_perf <- ggplot(perf_g, aes(y = performance, x = method, color = method)) +
  geom_boxplot(outlier.size = 0.4) + 
  scale_color_tron() + ylab(error_name) + xlab('') +
  theme(legend.position = 'None')

gg_perf
ggsave(filename = "simulation_study/figures/perf_fast.jpeg", plot = gg_perf, width = 3.5, height = 3)
