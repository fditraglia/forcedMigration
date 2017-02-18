solve_model <- function(V_cum_list, delta, tau_ell, tau_n, r, a0, a1,
                        land_parameters_list, gamma, alpha, n_cores = 1){

  solve_municipality_i <- function(i){
    p <- land_parameters_list[[i]]$p
    q <- land_parameters_list[[i]]$q
    H_bar <- land_parameters_list[[i]]$H_bar
    omega_n <- land_parameters_list[[i]]$omega_n
    V_cum <- V_cum_list[[i]]

    V_star_i <- get_V_star(delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                           omega_n, gamma, alpha)

    migration_cum_i <- get_migration_cum(V_cum, delta, tau_ell, tau_n, r, a0,
                                         a1, p, q, H_bar, omega_n)

    out <- list(V_star = V_star_i, migration_cum = migration_cum_i)
    return(out)
  }

  municipalities <- seq_along(land_parameters_list)
  model_solution <- parallel::mclapply(municipalities, solve_municipality_i,
                                       mc.cores = n_cores)

  V_total <- unlist(lapply(model_solution, function(x) x$V_star))
  D_cum <- lapply(model_solution, function(x) x$migration_cum)
  D_flow <- lapply(D_cum, function(x) diff(c(0, x)))

  out <- list(V_total = V_total, D_flow = D_flow)
  return(out)
}

solve_equilibrium_violence <- function(delta, tau_ell, tau_n, r, a0, a1,
                                       land_parameters_list, gamma, alpha,
                                       n_cores = 1){
  get_V_star_i <- function(i){
    p <- land_parameters_list[[i]]$p
    q <- land_parameters_list[[i]]$q
    H_bar <- land_parameters_list[[i]]$H_bar
    omega_n <- land_parameters_list[[i]]$omega_n
    get_V_star(delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n, gamma, alpha)
  }

  municipalities <- seq_along(land_parameters_list)
  V_total <- parallel::mclapply(municipalities, get_V_star_i, mc.cores = n_cores)
  return(unlist(V_total))
}


solve_equilibrium_migration_flow <- function(V_cum_list, delta, tau_ell, tau_n,
                                             r, a0, a1, land_parameters_list,
                                             n_cores = 1){
  get_migration_cum_i <- function(i){
    p <- land_parameters_list[[i]]$p
    q <- land_parameters_list[[i]]$q
    H_bar <- land_parameters_list[[i]]$H_bar
    omega_n <- land_parameters_list[[i]]$omega_n
    V_cum <- V_cum_list[[i]]
    get_migration_cum(V_cum, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n)
  }

  municipalities <- seq_along(land_parameters_list)
  D_cum <- parallel::mclapply(municipalities, get_migration_cum_i,
                              mc.cores = n_cores)

  D_flow <- lapply(D_cum, function(x) diff(c(0, x)))
  return(D_flow)
}


