demean_FE <- function(mat){
  t(apply(mat, 1, function(row) row - mean(row)))
}

bootstrap_F1 <- function(X, Y, n_bootstrap = 500,
                         q_lower = 0.05, q_upper = 0.95,
                         n_grid = 180){
  n_i <- nrow(X)
  n_T <- ncol(X)
  y_star <- c(t(demean_FE(Y)))
  X_star <- cbind(c(t(demean_FE(X))), c(t(demean_FE(X^2))))
  Q <- t(X)
  lower <- quantile(Q, q_lower)
  upper <- quantile(Q, q_upper)
  gamma_seq <- seq(lower, upper, length.out = n_grid)

  results_true <- get_S(X_star, y_star, Q, gamma_seq)
  S0 <- results_true$S0
  S1_gamma_hat <- results_true$S1_gamma_hat
  F1_true <- n_i * (n_T - 1) * (S0 - S1_gamma_hat) / S1_gamma_hat

  y_hat_0 <- X_star %*% results_true$coef0
  e_gamma_hat <- matrix(results_true$e_gamma_hat, n_T, n_i)

  bootstrap_draw <- function(){
    e_boot <- c(e_gamma_hat[, sample.int(n_i, replace = TRUE)])
    y_star_boot <- y_hat_0 + e_boot
    results_boot <- get_S(X_star, y_star_boot, Q, gamma_seq)
    S0_boot <- results_boot$S0
    S1_gamma_hat_boot <- results_boot$S1_gamma_hat
    return(n_i * (n_T - 1) * (S0_boot - S1_gamma_hat_boot) / S1_gamma_hat_boot)
  }

  F1_boot <- replicate(n_bootstrap, bootstrap_draw())
  out <- results_true
  out$F1_boot <- F1_boot
  out$F1 <- F1_true
  out$gamma_seq <- gamma_seq
  out$pvalue <- mean(F1_boot >= F1_true)
  out$n_i <- n_i
  out$n_T <- n_T
  return(out)
}

plot_gamma <- function(results, alpha = 0.05){
    c_alpha <- -2 * log(1 - sqrt(1 - alpha))
    S1_gamma <- results$S1
    S1_gamma_hat <- results$S1_gamma_hat
    gamma <- results$gamma_seq
    n_i <- results$n_i
    n_T <- results$n_T
    LR1 <- n_i * (n_T - 1) * (S1_gamma - S1_gamma_hat) / S1_gamma_hat
    plot(gamma, LR1, xlab = expression(gamma), ylab = expression(LR[1](gamma)),
         main = bquote(alpha == .(alpha)))
    abline(h = c_alpha, lty = 2, lwd = 2, col = 'red')
}

gamma_CI <- function(results, alpha = 0.05){
    c_alpha <- -2 * log(1 - sqrt(1 - alpha))
    S1_gamma <- results$S1
    S1_gamma_hat <- results$S1_gamma_hat
    gamma <- results$gamma_seq
    n_i <- results$n_i
    n_T <- results$n_T
    LR1 <- n_i * (n_T - 1) * (S1_gamma - S1_gamma_hat) / S1_gamma_hat
    return(gamma[which(LR1 < c_alpha)])
}

dgp_threshold <- function(d = 1, n_i = 200, n_t = 17, gamma = 6, b1 = 1,
                          b2 = -0.1){
  # Each row is an individual
  mu <- matrix(rnorm(n_i), n_i, n_t) # repeat these n_t times
  x <- matrix(runif(n_i * n_t, min = 0, max = 10), n_i, n_t)
  epsilon <- matrix(rnorm(n_i * n_t), n_i, n_t)
  z <- 1 * (x > gamma) # The true threshold indicator
  y <- mu + b1 * x + b2 * x^2 +  d * z + epsilon
  return(list(X = x, Y = y))
}
