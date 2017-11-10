# Solve the model for municipality i at parameters par
get_migration_flow_i <- function(i, par) {
  Vcum_i <- as.matrix(Vcum[i,])
  land_i <- land_parameters[i,]
  migration_cum <- get_migration_cum(Vcum_i,
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

get_dstar_friction <- function(par) {
  dstar_list <- lapply(1:nrow(Vcum), function(i) get_migration_flow_i(i, par))
  dstar <- do.call(rbind, dstar_list)
  colnames(dstar) <- colnames(Vcum)
  rownames(dstar) <- rownames(Vcum)
  dstar_lag <- cbind(rep(0, nrow(dstar)), dstar[,-ncol(dstar)])
  dstar_friction <- (1 - par$rho) * dstar + par$rho * dstar_lag
  return(dstar_friction)
}

negloglike_inner_D <- function(par_vec, Z, dstar_friction) {

  # Unpack parameters
  log_dbar <- par_vec[1]
  eta <- par_vec[2:17]
  nu <- par_vec[18:22]

  # Caclulate log-likelihood
  dbar <- exp(log_dbar)
  inner <- dbar + dstar_friction
  pop <- cross_section$popn1993
  mu <- (pop * t(t(inner) * exp(c(0, eta)))) %o% exp(nu)
  out <- -1 * (sum(Z * log(mu) - mu - lfactorial(Z)))

  # Calculate gradient
  dev <- Z - mu
  D_log_dbar <- dbar * sum(apply(dev, 3, function(mat) mat / inner))
  D_eta <- apply(dev, 2, sum)[-1] # The first time period doesn't have an eta
  D_nu <- apply(dev, 3, sum)
  attr(out, 'gradient') <- -1 * c(D_log_dbar, D_eta, D_nu)

  return(out)
}

negloglike_outer_D <- function(par_vec, Z) {

  par_model <- list(delta = par_vec[1],
                    tau_ell = par_vec[2],
                    tau_n = par_vec[3],
                    r = par_vec[4],
                    a0 = par_vec[5],
                    a1 = par_vec[6],
                    rho = par_vec[7])

  dstar_friction <- get_dstar_friction(par_model)
  fn <- function(x) negloglike_inner_D(x, Z, dstar_friction)
  suppressWarnings(opt <- nlm(f = fn, p = rep(0, 22)))
  return(opt$minimum)
}
