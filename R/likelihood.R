# Solve the model for municipality i at parameters par
get_migration_flow_i <- function(i, par) {
  Vcum_pop_i <- as.matrix(Vcum_pop[i,])
  land_i <- land_parameters[i,]
  migration_cum <- get_migration_cum(Vcum_pop_i,
                                     delta = par$delta,
                                     tau_ell = par$tau_ell,
                                     tau_n = par$tau_n,
                                     r = par$r,
                                     a0 = par$a0,
                                     a1 = par$a1,
                                     p = land_i$p,
                                     q = land_i$q,
                                     H_bar = land_i$H_bar,
                                     omega_n = land_i$omega_n)
  diff(c(0, migration_cum))
}

# Concentrated negative log-likelihood as a function of dbar and rho
negloglike_inner_D <- function(params, Z, dstar, dstar_lag) {
  dbar <- params[1]
  rho <- params[2]

  pop <- cross_section$popn1993
  D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
  sqrtZsum <- sqrt(sum(Z))
  Rt <- apply(Z, 2, sum) / sqrtZsum
  Rj <- apply(Z, 3, sum) / sqrtZsum
  Dsums <- colSums(D)
  mu <- t((t(D) / Dsums) * Rt) %o% Rj
  loglike <- (sum(Z * log(mu) - mu - lfactorial(Z)))

  return(-1 * loglike)
}

grad_inner_D <- function(params, Z, dstar, dstar_lag) {
  dbar <- params[1]
  rho <- params[2]

  # Calculate means \mu_{it}^j
  pop <- cross_section$popn1993
  D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
  sqrtZsum <- sqrt(sum(Z))
  Rt <- apply(Z, 2, sum) / sqrtZsum
  Rj <- apply(Z, 3, sum) / sqrtZsum
  Dsums <- colSums(D)
  mu <- t((t(D) / Dsums) * Rt) %o% Rj

  # Calculate gradient
  dev <- (Z - mu)
  onesJ <- rep(1, dim(Z)[3])
  Deriv_dbar <- sum((t(t(pop / D) - (sum(pop) / Dsums)) %o% onesJ * dev))
  popDiff <- pop * (dstar_lag - dstar)
  Deriv_rho <- sum(t(t(popDiff / D) - (colSums(popDiff) / Dsums)) %o% onesJ * dev)

  return(-1 * c(Deriv_dbar, Deriv_rho))
}

negloglike_outer_D <- function(par_vec, Z, return_inner = FALSE) {

  par_model <- list(delta = par_vec[1],
                    tau_ell = par_vec[2],
                    tau_n = par_vec[3],
                    r = par_vec[4],
                    a0 = par_vec[5],
                    a1 = par_vec[6])

  # Solve structural model
  f <- function(i) get_migration_flow_i(i, par_model)
  dstar <- do.call(rbind, lapply(1:nrow(Vcum_pop), f))
  dstar_lag <- cbind(rep(0, nrow(dstar)), dstar[,-ncol(dstar)])

  # Maximize concentrated log-likelihood over dbar and rho
  f <- function(x) negloglike_inner_D(x, Z, dstar, dstar_lag)
  Deriv_f <- function(x) grad_inner_D(x, Z, dstar, dstar_lag)
  opt <- optim(c(0.01, 0.1), fn = f,  gr = Deriv_f,
               lower = c(1e-10, 0), upper = c(1, 1), method = 'L-BFGS-B')

  # Optionally return the optimal values of the "inner" parameters dbar, rho,
  # eta, and nu. (The default is *not* to return these.)
  if(return_inner) {
    dbar <- opt$par[1]
    rho <- opt$par[2]
    negloglike <- opt$value

    # Calculate the optimal values of eta and nu
    pop <- cross_section$popn1993
    D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
    log_S <- log(colSums(D))
    log_F <- log(apply(Z, 2, sum)) -  log(sum(Z))
    eta <- log_F - log_F[1] + log_S[1] - log_S
    nu <- log(apply(Z, 3, sum)) + log_F[1] - log_S[1]
    return(list(negloglike = negloglike, dbar = dbar, rho = rho,
                eta = eta[-1], nu = nu)) # First element of eta normalized to 0
  } else {
    return(opt$value)
  }
}

# Payoffs for contract as a function of migration parameters
get_payoffs_list <- function(delta, tau_ell, tau_n, r, a0, a1,
                             n_cores = 1) {
  get_payoffs_i <- function(i) {
    get_payoffs(delta, tau_ell, tau_n, r, a0, a1,
                land_parameters[i,]$p,
                land_parameters[i,]$q,
                land_parameters[i,]$H_bar,
                land_parameters[i,]$omega_n,
                cross_section[i,]$popn1993)
  }
  parallel::mclapply(1:nrow(land_parameters), get_payoffs_i, mc.cores = n_cores)
}

# Contracts: lambda_star and S_e_star for each municipality as a function of payoffs
# and contract parameters
get_contracts <- function(gamma, alpha, payoffs) {
  get_contract_i <- function(i) {
    get_contract(gamma, alpha,
                  payoffs[[i]]$D_max,
                  payoffs[[i]]$X_max,
                  payoffs[[i]]$Dstar,
                  payoffs[[i]]$Xstar,
                  cross_section$num_families[i])
  }
  contracts <- lapply(seq_along(payoffs), get_contract_i)
  lambda_star <- sapply(contracts, function(x) x$lambda_star)
  S_e_star <- sapply(contracts, function(x) x$S_e_star)
  upper <- sapply(contracts, function(x) x$lambda_upper)
  out <- data.frame(lam = lambda_star, S_e = S_e_star, upper = upper)
  return(out)
}

neg_loglike_inner_V <- function(params, V, covariates, contracts) {
  # Note: the covariates are assumed *not* to contain an intercept
  #       so we add one below.

  lam <- contracts$lam
  S_e <- contracts$S_e
  # We define y = S_e - x'beta, hence the negative signs in X
  X <- cbind(S_e,  rep(-1, length(S_e)), -1 * covariates)
  y <- X %*% params

  Vzero <- V == 0
  Vpos <- V != 0
  out1 <- -1 * sum(log(1 + exp(y[Vzero])))
  out2 <- -1 * sum(log(1 + exp(-y[Vpos])) + log(exp(lam[Vpos]) - 1) -
                      V_obs[Vpos] * log(lam[Vpos]))
  out3 <- -1 * sum(log(factorial(V_obs[Vpos]))) # This term is a constant
  out <- -1 * (out1 + out2 + out3)

  grad_zero <- -1 * t(X[Vzero,]) %*% (1 / (1 + exp(-y[Vzero])))
  grad_pos <- t(X[Vpos,]) %*%  (1 / (1 + exp(y[Vpos])))
  attr(out, 'gradient') <- -1 * (grad_zero + grad_pos)

  return(out)
}

# The argument covariates (a matrix) is assumed to *exclude* a constant term.
# The constant is added by neg_loglike_inner_V
neg_loglike_outer_V <- function(par_vec, V, covariates, payoffs) {
  par_model <- list(gamma = par_vec[1],
                    alpha = par_vec[2])

  contracts <- get_contracts(gamma = par_model$gamma,
                             alpha = par_model$alpha, payoffs)
  f <- function(x) neg_loglike_inner_V(x, V, covariates, contracts)
  opt <- suppressWarnings(nlm(f, p = rep(0, ncol(covariates) + 2),
                              check.analyticals = FALSE))
  return(opt$minimum)
}
