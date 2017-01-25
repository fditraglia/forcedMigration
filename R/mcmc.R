moment_condition <- function(delta, tau_ell, tau_n, r, a0, a1,
                             land_parameters_list, gamma, alpha,
                             V_total_obs, V_cum_obs, D_flow_obs,
                             n_cores = 1){
  # Solve model
  V_total_model <- solve_equilibrium_violence(delta, tau_ell, tau_n, r, a0, a1,
                                              land_parameters_list, gamma, alpha,
                                              n_cores)
  D_flow_model <- solve_equilibrium_migration_flow(V_cum_obs, delta, tau_ell,
                                                   tau_n, r, a0, a1,
                                                   land_parameters_list,
                                                   n_cores)
  # Convert migration flows to matrix
  n_periods <- length(D_flow_obs[[1]])
  matD_flow_obs <- matrix(unlist(D_flow_obs), ncol = n_periods, byrow = TRUE)
  matD_flow_model <- matrix(unlist(D_flow_model), ncol = n_periods, byrow = TRUE)

  # Construct moment conditions
  D_MC <- colMeans(matD_flow_model - matD_flow_obs)
  V_MC <- mean(V_total_model - V_total_obs)
  MC <- c(D_MC, V_MC)
  return(sum(MC^2))
}

# Hard-code the simulation data into this moment condition to save typing
moment_condition_sim <- function(delta, tau_ell, tau_n, r, a0, a1, gamma,
                                 alpha, n_cores = 1){

  out <- moment_condition(delta = delta, tau_ell = tau_ell,
                          tau_n = tau_n, r = r, a0 = a0, a1 = a1,
                          land_parameters_list = sims$land_params,
                          gamma = gamma, alpha = alpha,
                          V_total_obs = sims$V_total,
                          V_cum_obs = sims$V_cum,
                          D_flow_obs = sims$D_flow_obs,
                          n_cores = n_cores)
  return(out)
}


draw_CH_chain_simdata <- function(n_draws, n_cores){

  true_params <- unlist(sims$model_params)
  lower <- true_params / c(10, 1, 1, 1, 1, 1, 10, 10)
  upper <- true_params * c(10, 1, 1, 1, 1, 1, 10, 10)
  starting_values <- true_params + c(0.8, 0, 0, 0, 0, 0, 0.75, 0.005)
  sigma <- (upper - lower) / 3
  n_params <- length(sigma)

  draws <- matrix(NA_real_, nrow = n_draws + 1, ncol = n_params)
  colnames(draws) <- names(sims$model_params)
  Like <- rep(NA, n_draws)

  accept_count <- 0L
  draws[1,] <- starting_values
  Like[1] <- -1 * do.call(moment_condition_sim, c(as.list(starting_values),
                                                  list(n_cores = n_cores)))
  for(i in 2:(n_draws + 1)){
    proposal <- draws[i-1,] + rnorm(n_params, mean = 0, sd = sigma)

    if(any(proposal > upper) || any(proposal < lower)){
      draws[i,] <- draws[i-1,]
      Like[i] <- Like[i-1]
    } else {
      Like_proposal <- -1 * do.call(moment_condition_sim, c(as.list(proposal),
                                                      list(n_cores = n_cores)))
      rho <- min(exp(Like_proposal - Like[i-1]), 1)
      U <- runif(1)
      if(U < rho){
        draws[i,] <- proposal
        Like[i] <- Like_proposal
        accept_count <- accept_count + 1
      } else {
        draws[i,] <- draws[i-1,]
        Like[i] <- Like[i-1]
      }
    }
  }
  call <- list(start = starting_values, sigma = sigma, lower = lower,
               uppwer = upper)
  out <- list(draws = draws, Like = Like, accept_rate = accept_count / n_draws,
              call = call)
  return(out)
}



