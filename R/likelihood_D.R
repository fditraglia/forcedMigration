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
negloglike_inner_D <- function(params, Z, dstar, dstar_lag, time_effects = TRUE) {
  dbar <- params[1]
  rho <- params[2]

  pop <- cross_section$popn1993
  D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
  if(time_effects) {
    sqrtZsum <- sqrt(sum(Z))
    Rt <- apply(Z, 2, sum) / sqrtZsum
    Rj <- apply(Z, 3, sum) / sqrtZsum
    Dsums <- colSums(D)
    mu <- t((t(D) / Dsums) * Rt) %o% Rj
  } else {
    mu <- (D / sum(D)) %o% apply(Z, 3, sum)
  }
  loglike <- (sum(Z * log(mu) - mu - lfactorial(Z)))

  return(-1 * loglike)
}

grad_inner_D <- function(params, Z, dstar, dstar_lag, time_effects = TRUE) {
  dbar <- params[1]
  rho <- params[2]

  pop <- cross_section$popn1993
  D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
  popDiff <- pop * (dstar_lag - dstar)
  onesJ <- rep(1, dim(Z)[3])

  if(time_effects) {

    sqrtZsum <- sqrt(sum(Z))
    Rt <- apply(Z, 2, sum) / sqrtZsum
    Rj <- apply(Z, 3, sum) / sqrtZsum
    Dsums <- colSums(D)
    mu <- t((t(D) / Dsums) * Rt) %o% Rj
    dev <- (Z - mu)
    Deriv_dbar <- sum(t(t(pop / D) - (sum(pop) / Dsums)) %o% onesJ * dev)
    Deriv_rho <- sum(t(t(popDiff / D) - (colSums(popDiff) / Dsums)) %o% onesJ * dev)

  } else {

    Dsum <- sum(D)
    mu <- (D / Dsum) %o% apply(Z, 3, sum)
    dev <- (Z - mu)
    Deriv_dbar <- sum((pop / D - dim(D)[2] * sum(pop) / Dsum) %o% onesJ * dev)
    Deriv_rho <- sum((popDiff / D - sum(popDiff) / Dsum) %o% onesJ * dev)

  }
  return(-1 * c(Deriv_dbar, Deriv_rho))
}

negloglike_outer_D <- function(par_vec, Z, X = NULL, time_effects = TRUE,
                               return_inner = FALSE, ncores = 1) {

  # Case of no covariates (X): identical to what we had before so old code still works
  # NOTE: in this case the parameters are *not* within an exp!
  if(is.null(X)) {
    stopifnot(length(par_vec) == 6)
    par_model <- list(delta = par_vec[1],
                      tau_ell = par_vec[2],
                      tau_n = par_vec[3],
                      r = par_vec[4],
                      a0 = par_vec[5],
                      a1 = par_vec[6])

    # Solve structural model with common parameters for each municipality
    solve_model_i <- function(i) get_migration_flow_i(i, par_model)

  } else {
    # In this case, the parameters are within an exp! See the comments in
    # likelihood_helper.R for more information

    stopifnot(all.equal(names(X), c('delta', 'tau_ell', 'tau_n', 'r', 'a0', 'a1')))
    heterog_pars <- get_heterog_pars(par_vec, X)

    # Function to solve structural model
    solve_model_i <- function(i) {
      par_model_i <- list(delta = heterog_pars$delta[i],
                          tau_ell = heterog_pars$tau_ell[i],
                          tau_n = heterog_pars$tau_n[i],
                          r = heterog_pars$r[i],
                          a0 = heterog_pars$a0[i],
                          a1 = heterog_pars$a1[i])
      get_migration_flow_i(i, par_model_i)
    }
  } # END else

  dstar <- do.call(rbind, parallel::mclapply(1:nrow(Vcum_pop),
                                             solve_model_i, mc.cores = ncores))
  dstar_lag <- cbind(rep(0, nrow(dstar)), dstar[,-ncol(dstar)])

  # Maximize concentrated log-likelihood over dbar and rho
  f <- function(x) negloglike_inner_D(x, Z, dstar, dstar_lag,
                                      time_effects = time_effects)
  Deriv_f <- function(x) grad_inner_D(x, Z, dstar, dstar_lag,
                                      time_effects = time_effects)
  opt <- optim(c(0.01, 0.1), fn = f,  gr = Deriv_f,
               lower = c(1e-10, 0), upper = c(1, 1), method = 'L-BFGS-B')

  # Optionally return the optimal values of the "inner" parameters dbar, rho,
  # eta, and nu. (The default is *not* to return these.) Also return the
  # heterogenous parameters
  if(return_inner) {
    dbar <- opt$par[1]
    rho <- opt$par[2]
    negloglike <- opt$value

    # Calculate the optimal values of eta and nu
    pop <- cross_section$popn1993
    D <- pop * (dbar + (1 - rho) * dstar + rho * dstar_lag)
    if(time_effects) {
      log_S <- log(colSums(D))
      log_F <- log(apply(Z, 2, sum)) -  log(sum(Z))
      eta <- log_F - log_F[1] + log_S[1] - log_S
      nu <- log(apply(Z, 3, sum)) + log_F[1] - log_S[1]
    } else {
      eta <- rep(0, dim(Z)[2])
      nu <- log(apply(Z, 3, sum)) - log(sum(D))
    }

    out <- list(negloglike = negloglike, dbar = dbar, rho = rho,
                eta = eta[-1], nu = nu, D = D, dstar = dstar,
                dstar_lag = dstar_lag, pop = pop)

    if(!is.null(X)) {
      out$X <- X
      out$parvecs <- heterog_pars
    }
    return(out)

  } else {
    return(opt$value)
  }
}
