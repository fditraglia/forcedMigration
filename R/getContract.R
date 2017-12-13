get_contractR <- function(gamma, alpha, D_max, X_max, Dstar, Xstar, nfam) {
  # Let B = nfam * (X - gamma * D), i.e. surplus plus the cost of violence
  B_max <- nfam * (X_max - gamma * D_max)
  Bstar <- nfam * (Xstar - gamma * Dstar)

  # Expected surplus as a function of lambda
  S_e <- function(lambda) {
    V_seq <- seq_along(Bstar)
    denom <- 1 - exp(-lambda)
    probs <- dpois(V_seq, lambda) / denom
    sum(probs * Bstar) + (1 - sum(probs)) * B_max -
      0.5 * alpha * lambda * (1 + lambda) / denom
  }
  # S_e <- function(lambda) {
  # V_seq <- seq_along(Bstar)
  # log_denom <- log(1.0 - exp(-lambda))
  # probs <- exp(dpois(V_seq, lambda, log = TRUE) - log_denom)
  # ExpectedSurplus <- sum(probs * Bstar) + (1.0 - sum(probs)) * B_max
  #     - 0.5 * alpha * exp(log(lambda) + log(1.0 + lambda) - log_denom)
  #  return(ExpectedSurplus)
  # }

  # Maximize expected surplus using the Brent algorithm
  lower <- 0.1; # Smallest lambda allowed: corresponds to approx. 95%
                      # probability of drawing a 1 from a zero-truncated
                      # Poisson distribution
  S_e_lower <- S_e(lower) # Expected surplus at lamba = lower
  M <- max(max(Bstar), B_max)
  upper <- 2.0 * (M - S_e_lower) / alpha
  opt <- optimize(f = S_e, interval = c(lower, upper), maximum = TRUE)
  lambda_star <- opt$maximum
  S_e_star <- opt$objective

  # Given the shape of the objective function, the only way for Brent to fail
  # is if it finds a local optimum that gives a lower expected surplus than
  # setting lambda equal to lower.
  if(S_e_star < S_e_lower) {
    lambda_star = lower
    S_e_star = S_e_lower
  }

  return(list(lambda_upper = upper,
              lambda_star = lambda_star,
              S_e_star = S_e_star))
}
