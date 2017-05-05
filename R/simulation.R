sticky_bins <- function(n_balls, n_bins, k = 1){
  stopifnot(is.integer(n_balls) && is.integer(n_bins))
  stopifnot(all.equal(length(n_balls), 1L) && all.equal(length(n_bins), 1L))
  stopifnot((n_balls >= 0L) && (n_bins > 0L))
  bins <- rep(0L, n_bins)
  for(i in seq_len(n_balls)){
    bin_probs <- (bins + 1L) / length(bins)
    next_bin <- sample(seq_len(n_bins), size = 1, prob = k * bins + 1)
    bins[next_bin] <- bins[next_bin] + 1L
  }
  return(bins)
}

simulate_from_model <- function(delta, tau_ell, tau_n, r, a0, a1,
                                land_parameters_list, gamma, alpha,
                                num_families,
                                d_intercept = 0.0003,
                                # d_intercept gives Poisson mean of about 0.7
                                # when scaled up by median num_families (2330)
                                V_intercept = 0.2,
                                # V_intercept is already on the count scale
                                # since violence is not per-capita
                                n_cores = 1){

  V_total <- solve_equilibrium_violence(delta, tau_ell, tau_n, r, a0, a1,
                                        land_parameters_list, gamma, alpha,
                                        n_cores)

  V_total_obs <- rpois(length(V_total), V_intercept + V_total)
  chop_violence <- function(x) cumsum(sticky_bins(x, 17L, 3))
  V_cum_list <- lapply(V_total_obs, chop_violence)
  D_flow <- solve_equilibrium_migration_flow(V_cum_list, delta, tau_ell, tau_n,
                                             r, a0, a1, land_parameters_list,
                                             n_cores)
  d_noise <- function(d, n){
    rpois(length(d), n * (d + d_intercept)) / n
  }
  D_flow_obs <- Map(d_noise, D_flow, num_families)

  out <- list(V_total = V_total, V_total_obs = V_total_obs,
              V_cum_obs = V_cum_list,
              D_flow = D_flow, D_flow_obs = D_flow_obs)
}

