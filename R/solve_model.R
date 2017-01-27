solve_equilibrium_violence <- function(delta, tau_ell, tau_n, r, a0, a1,
                                       land_parameters_list, gamma, alpha,
                                       cluster = NULL, n_cores = 1){

  migration_params <- list(delta = delta, tau_ell = tau_ell, tau_n = tau_n,
                           r = r, a0 = a0, a1 = a1)
  surplus_params <- list(gamma = gamma, alpha = alpha)

  get_V_star_i <- function(i){
    params <- c(migration_params, surplus_params, land_parameters_list[[i]])
    do.call(get_V_star, params)
  }

  if(!is.null(cluster)){
    V_total <- parallel::clusterApplyLB(cluster, seq_along(land_parameters_list),
                                   get_V_star_i)
  } else {
    V_total <- parallel::mclapply(seq_along(land_parameters_list),
                                  get_V_star_i, mc.cores = n_cores)
  }
  return(unlist(V_total))
}


solve_equilibrium_migration_flow <- function(V_cum_list, delta, tau_ell, tau_n,
                                             r, a0, a1, land_parameters_list,
                                             cluster = NULL, n_cores = 1){

  migration_params <- list(delta = delta, tau_ell = tau_ell, tau_n = tau_n,
                           r = r, a0 = a0, a1 = a1)
  get_migration_cum_i <- function(i){
    V_list <- list(V_cum = V_cum_list[[i]])
    params <- c(V_list, migration_params, land_parameters_list[[i]])
    do.call(get_migration_cum, params)
  }

  if(!is.null(cluster)){
    D_cum <- parallel::clusterApplyLB(cluster, seq_along(land_parameters_list),
                                 get_migration_cum_i)
  } else {
    D_cum <- parallel::mclapply(seq_along(land_parameters_list),
                              get_migration_cum_i, mc.cores = n_cores)
  }
  D_flow <- lapply(D_cum, function(x) diff(c(0, x)))
  return(D_flow)
}


