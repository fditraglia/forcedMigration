# Payoffs for contract as a function of migration parameters
get_payoffs_list <- function(par_vec, X = NULL, n_cores = 1) {

  # Case of no covariates (X) in the migration problem
  # NOTE: in this case the parameters are *not* within an exp!
  if (is.null(X)) {

    stopifnot(length(parvec) == 6)
    delta <- par_vec[1]
    tau_ell <- par_vec[2]
    tau_n <- par_vec[3]
    r <- par_vec[4]
    a0 <- par_vec[5]
    a1 <- par_vec[6]

    # Function to calculate payoffs for a given municipality
    get_payoffs_i <- function(i) {
      get_payoffs(delta, tau_ell, tau_n, r, a0, a1,
                  land_parameters[i,]$p,
                  land_parameters[i,]$q,
                  land_parameters[i,]$H_bar,
                  land_parameters[i,]$omega_n,
                  cross_section[i,]$popn1993)
    }

  } else {
    # NOTE: in this case the parameters are *are* within an exp!
    stopifnot(all.equal(names(X), c('delta', 'tau_ell', 'tau_n', 'r', 'a0', 'a1')))
    heterog_pars <- get_heterog_pars(par_vec, X)

    # Function to calculate payoffs for a given municipality
    get_payoffs_i <- function(i) {
      get_payoffs(delta = heterog_pars$delta[i],
                  tau_ell = heterog_pars$tau_ell[i],
                  tau_n = heterog_pars$tau_n[i],
                  r = heterog_pars$r[i],
                  a0 = heterog_pars$a0[i],
                  a1 = heterog_pars$a1[i],
                  land_parameters[i,]$p,
                  land_parameters[i,]$q,
                  land_parameters[i,]$H_bar,
                  land_parameters[i,]$omega_n,
                  cross_section[i,]$popn1993)
    }
  }
  parallel::mclapply(1:nrow(land_parameters), get_payoffs_i, mc.cores = n_cores)
}


#-------------------------------------------------------------------------------
# NOTE 2018-06-14: this needs to be changed so that alpha and gamma can depend
# on covariates. Other than that, it can stay the same.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# NOTE 2018-06-14: this is the "inner" logit regression, nothing needs to change!
# But we may want to rename the argument "covariates" to be clearer about what
# covariates we mean exactly, e.g. "inner_V_covariates" or cost_covariates or
# something like that.
#-------------------------------------------------------------------------------
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

#-------------------------------------------------------------------------------
# NOTE 2018-06-14: this only needs to change in terms of what arguments it takes
# and how it unpacks them. Need to think carefully about the names of different
# things here, e.g. "covariates"
#-------------------------------------------------------------------------------
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
