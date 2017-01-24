chop_violence <- function(V_total, n_periods = 17){
  stopifnot(is.atomic(V_total))
  n_sample <- length(V_total)
  weights <- rep(1, n_periods)
  dirichlet_shape <- rep(weights, times = n_sample)
  gamma_draws <- rgamma(n_periods * n_sample, shape = dirichlet_shape, scale = 1)
  gamma_draws <- matrix(gamma_draws, nrow = n_sample, ncol = n_periods)
  dirichlet_draws <- t(apply(gamma_draws, 1, function(x) x / sum(x)))
  V_seq <- V_total * dirichlet_draws
  V_cum <- t(apply(V_seq, 1, cumsum))
  V_cum_list <- as.list(data.frame(t(V_cum)))
  names(V_cum_list) <- NULL
  return(V_cum_list)
}

simulate_from_model <- function(delta, tau_ell, tau_n, r, a0, a1,
                                land_parameters_list, gamma, alpha,
                                n_cores = 1){

  V_total <- solve_equilibrium_violence(delta, tau_ell, tau_n, r, a0, a1,
                                        land_parameters_list, gamma, alpha,
                                        n_cores = n_cores)
  V_cum_list <- chop_violence(V_total)
  D_flow <- solve_equilibrium_migration_flow(V_cum_list, delta, tau_ell, tau_n,
                                             r, a0, a1, land_parameters_list,
                                             n_cores = n_cores)

  add_noise <- function(x) x + rnorm(length(x), mean = 0, sd = x / 3)
  V_total_obs <- add_noise(V_total)
  D_flow_obs <- lapply(D_flow, add_noise)
  out <- list(V_total = V_total, V_total_obs = V_total_obs,
              D_flow = D_flow, D_flow_obs = D_flow_obs)
}

