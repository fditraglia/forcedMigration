demean_FE <- function(mat){
  t(apply(mat, 1, function(row) row - mean(row)))
}

#' Bootstrap test for threshold effect following Hansen (1999, JOE)
#'
#' @param X Panel data for x. Assumes that each row is an individual and each
#' column is a time period.
#' @param Y Panel data for y. Assumes that each row is an individual and each
#' column is a time period.
#' @param n_bootstrap Number of bootstrap draws
#' @param q_lower Quantile of x above which we will search for the threshold
#' @param q_upper Quantile of of x below which we will search for the threshold
#' @param n_grid  Number of grid points between \code{q_lower} and \code{q_lower}
#' over which we will search for the the threshold
#'
#' @return List with bootstrap results:
#' \itemize{
#'    \item \code{gamma_hat} Estimated threshold location
#'    \item \code{gamma_seq} Grid of values considered for threshold location
#'    \item \code{S0} Sum of squared residuals for restricted model with no threshold
#'    \item \code{coef0} Estimated regression coefficients for restricted model with no
#'    \item \code{e_tilde} Residuals from the restricted model with no threshold
#'    threshold
#'    \item \code{e_gamma_hat} Residuals for threshold model with gamma set equal to
#'    gamma_hat
#'    \item \code{S1_gamma_hat} Sum of squared residuals for threshold model with gamma
#'    set equal to gamma_hat
#'    \item \code{S1} Vector corresponding to gamma_seq. Each element is the sum of squared
#'    residuals for a thresdhold model with that particular value of gamma.
#'    \item \code{F1} Value of the F-test statistic for the input dataset
#'    \item \code{F1_boot} Vector of values for F-test statistic in each of the
#'    bootstrap samples
#'    \item \code{pvalue} Bootstrap p-value for test of the null hypothesis that there is no
#'    threshold, calculated from F1 and F1_boot
#'    \item \code{n_i} Number of individuals in the panel
#'    \item \code{n_T} Number of time periods in the panel
#'  }
#'
#' @examples
#' test_data <- dgp_threshold()
#' results <- bootstrap_F1(test_data$X, test_data$Y)
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

#' Plot confidence set for estimated threshold location following Hansen (1999, JOE)
#'
#' @param results List of results from bootstrap_F1
#' @param alpha Confidence level. Defaults to 0.05
#'
#' @return Plots the results of the likelihood ratio test from \code{gamma_CI}.
#' Each point is a test statistic corresponding to a particular value of
#' \code{results$gamma_seq}. The red line gives the critical value for a particular
#' choice of \code{alpha}, which defaults to 0.05. Points below the red line
#' correspond to values of gamma that cannot be rejected.
#'
#' @examples
#' test_data <- dgp_threshold()
#' results <- bootstrap_F1(test_data$X, test_data$Y)
#' plot_gamma(results)
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

#' Confidence set for estimated threshold location following Hansen (1999, JOE)
#'
#' @param results List of results from bootstrap_F1
#' @param alpha Confidence level. Defaults to 0.05
#'
#' @return Vector of values for gamma from \code{results$gamma_seq} that do not
#' give a significantly different sum of squared residuals than \code{results$gamma_hat}
#' based on a likelihood ratio test with level alpha
#'
#' @examples
#' test_data <- dgp_threshold()
#' results <- bootstrap_F1(test_data$X, test_data$Y)
#' gamma_CI(results)
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

#' Generate test dataset to test for threshold effects
#'
#' @param d Size of the threshold effect
#' @param n_i Number of individuals
#' @param n_t Number of time periods
#' @param gamma Location of threshold
#' @param b1 Coefficient on x
#' @param b2 Coefficient on x-squared
#'
#' @return List with panel data X and Y formatted appropriately for bootstrap_F1
#'
#' @examples
#' test_data <- dgp_threshold()
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
