library(EQL)
library(SDForest)
library(ggplot2)
library(gridExtra)
library(ggsci)
library(tidyr)

mc.cores <- 1

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

f_hermite <- function(x, beta, js){
    # function to generate f_X
    # x: covariates
    # beta: parameter vector
    # js: relevant covariates

    # number of relevant covariates
    m <- length(js)

    # complexity of f_X
    complexity <- length(beta) / m

    if(is.null(dim(beta))) beta <- matrix(beta)

    # calculate f_X
    do.call(sum, lapply(1:m, function(i) {
        j <- js[i]
        # select beta for covariate j
        res <- lapply(1:complexity, function(k) {
            beta[k, i] * hermite(x[j], k-1)
        })
        Reduce('+', res)
    }))
}

simulate_nonlinear_confounding <- function(q, p, n, m, df){
  H <- matrix(rnorm(n*q), ncol=q)
  
  Betas <- replicate(p, matrix(runif(q*df, -1, 1), nrow = df), simplify = F)
  X <- sapply(Betas, function(beta) {apply(H, 1, function(h) f_hermite(h, beta, 1:q))})
  X <- X + matrix(rnorm(n*p), ncol = p)

  complexity <- 2
  # random parameter for fourier basis
  beta <- runif(m * complexity * 2, -1, 1)
  # random sparse subset of covariates in X
  js <- sample(1:p, m)

  # random coefficient vector delta
  delta <- rnorm(q, 0, 2)

  # generate f_X
  f_X <- apply(X, 1, function(x) f_four(x, beta, js))
  Y <- f_X + H %*% delta + rnorm(n, 0, 0.01)
  return(list(X = X, Y = Y, f_X = f_X, j = js))
}

# Visualization of hermite-polynomials
p <- 1
q <- 1
n <- 1000
max_df <- 4
  
set.seed(42)
H <- matrix(rnorm(n*q), ncol=q)
X <- sapply(1:max_df, function(df){
  Betas <- replicate(p, matrix(runif(q*df, -1, 1), nrow = df), simplify = F)
  X <- sapply(Betas, function(beta) {apply(H, 1, function(h) f_hermite(h, beta, 1:q))})
})

X <- as.vector(X)
H <- rep(H, max_df)
df <- as.factor(rep(paste('df =', 1:max_df), each = n))

herm_df <- data.frame(X, H, df)
gg_herm <- ggplot(herm_df, aes(x = H, y = X)) + 
  geom_line() + 
  facet_grid(~df) + 
  theme_bw()
gg_herm
ggsave(filename = "simulation_study/figures_nl/herm.jpeg", plot = gg_herm, width = 10, height = 4)

# Example simulations
set.seed(42)
p <- 30
q <- 1
n <- 50
df <- 4
m <- 1

dat <- simulate_nonlinear_confounding(q = q, p = p, n = n, m = m, df = df)
X <- dat$X
Y <- dat$Y
f_X <- dat$f_X
js <- dat$j

# spiking of df singular values
Q <- get_Q(X, 'trim')
d <- svd(X)$d
sing <- data.frame(d, i = 1:p)

# transformed response and f(X)
df <- data.frame(f = f_X, Y, QY = Q %*% Y, Qf = Q %*% f_X)

# Estimation of models
fit <- SDForest(x = X, y = Y, mc.cores = mc.cores)
dep <- partDependence(fit, js[1], mc.cores = mc.cores)
fit2 <- SDForest(x = X, y = Y, Q_type = 'no_deconfounding', mc.cores = mc.cores)
dep2 <- partDependence(fit2, js[1], mc.cores = mc.cores)

# Comparison of Variable importance
imp_1 <- fit$var_importance / max(fit$var_importance)
imp_2 <- fit2$var_importance / max(fit2$var_importance)
true_imp <- rep('spurious     ', length(imp_1))
true_imp[js] <- 'causal'

imp_data <- data.frame(SDF = imp_1, ranger = imp_2, Covariates = as.factor(true_imp))

# Quantitative performance comparison
error_name <- expression("||"*f^0*(x[test]) - widehat(f(x[test]))*"||"[2]^2 / n[test])

performance_measure <- function(n, p, q, n_test){
  data <- simulate_nonlinear_confounding(q = q, p = p, n = n+n_test, m = 4, df = 4)
  data_train <- data.frame(data$X[1:n,], Y = data$Y[1:n])
  data_test <- data.frame(data$X[(n+1):(n+n_test),], Y = data$Y[(n+1):(n+n_test)])
  
  fit <- SDForest(Y ~ ., data_train, mc.cores = mc.cores)
  fit2 <- SDForest(Y ~ ., data_train, Q_type = "no_deconfounding", mc.cores = mc.cores)
  
  pred <- predict(fit, data_test)
  pred2 <- predict(fit2, data_test)
  
  mse <- (data$f_X[(n+1):(n+n_test)] - pred)
  mse2 <- (data$f_X[(n+1):(n+n_test)] - pred2)

  return(list(SDF = mse, RF= mse2))
}

perf <- replicate(10, sapply(performance_measure(n, p, q, 500), 
                              function(x)mean(x**2)))
perf <- t(perf)

perf_g <- gather(data.frame(perf), key = 'method', value = 'performance')
save(perf_g, fit, fit2, dep, dep2,
     file = "simulation_study/results/nonlin_confounding.RData")

# Plots
gg_sing <- ggplot(sing, aes(x = i, y = d)) + 
  geom_point(aes()) + ylab(expression(lambda[i])) + 
  theme_bw()
gg_sing
ggsave(filename = "simulation_study/figures_nl/sing.jpeg", plot = gg_sing, width = 7, height = 3)

plain <- ggplot(df, aes(y = f, x = Y)) + 
  geom_point(size = 0.4) + theme_bw() + 
  ylab("f(X)")

transformed <- ggplot(df, aes(y = Qf, x = QY)) + 
  geom_point(size = 0.4) + theme_bw() + 
  ylab("Qf(X)")

pt <- grid.arrange(plain, transformed)
ggsave(filename = "simulation_study/figures_nl/pt_nl.jpeg", plot = pt, width = 5, height = 10)

ggimp <- ggplot(imp_data, aes(x = SDF, y = ranger, col = Covariates)) + 
  geom_point(size = 0.5) + theme_bw() + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Normalized to [0, 1]') + 
  theme(legend.title = element_blank())

ggimp_log <- ggplot(imp_data, aes(x = log(SDF), y = log(ranger), col = Covariates)) + 
  geom_point(size = 0.5) + theme_bw() + xlab('') + 
  ylab('') + scale_color_tron() + ggtitle('Logarithmic scale') + 
  theme(legend.title = element_blank())

gg_imp <- grid_arrange_shared_legend(ggimp, ggimp_log, position = 'right',
                                     left = 'Variable importance ranger', 
                                     bottom = 'Variable importance SDForest')
ggsave(filename = "simulation_study/figures_nl/imp_nl.jpeg", plot = gg_imp, width = 10, height = 4)

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
                                         rf = "no deconfounding"))
gg_dep
ggsave(filename = "simulation_study/figures_nl/dep_nl.jpeg", plot = gg_dep, width = 5, height = 5)


gg_perf <- ggplot(perf_g, aes(y = performance, x = method, fill = method)) +
  geom_boxplot(outlier.size = 0.4) + theme_bw() + 
  scale_fill_tron() + ylab(error_name) + xlab('') +
  theme(legend.position = 'None')

gg_perf
ggsave(filename = "simulation_study/figures_nl/perf_nl.jpeg", plot = gg_perf, width = 5, height = 5)
