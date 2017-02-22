moment_condition <- function(delta, tau_ell, tau_n, r, a0, a1,
                             land_parameters_list, gamma, alpha,
                             V_total_obs, V_cum_obs, D_flow_obs,
                             n_cores = 1){
  # Solves the model at specified parameter values and constructs matrix MC each
  # column of which is an m_i(theta) for some municipality i
  model_solution <- solve_model(V_cum_list = V_cum_obs, delta,
                                tau_ell, tau_n, r, a0, a1, land_parameters_list,
                                gamma, alpha, n_cores)
  V_total_model <- model_solution$V_total
  D_flow_model <- model_solution$D_flow

  # Convert migration flows to matrix
  n_periods <- length(D_flow_obs[[1]])
  matD_flow_obs <- matrix(unlist(D_flow_obs), ncol = n_periods, byrow = TRUE)
  matD_flow_model <- matrix(unlist(D_flow_model), ncol = n_periods, byrow = TRUE)

  # Construct moment functions m_i(theta) for each municipality
  D_MC <- matD_flow_model - matD_flow_obs
  colnames(D_MC) <- paste0('D_MC_t', 1:ncol(D_MC))
  V_MC <- V_total_model - V_total_obs
  MC <- t(cbind(V_MC, D_MC))
  return(MC)
}

Like_CUE <- function(MC){
  # Evaluates the pseudo-likelihood based on the continuous-updating GMM
  # criterion function, where each row of MC is an m_i(theta) for some
  # municipality i
  n <- ncol(MC)
  g <- rowMeans(MC)
  Omega <- var(t(MC))
  L <- t(chol(Omega))
  g_weighted <- forwardsolve(L, g)
  return(-0.5 * n * sum(g_weighted^2))
}

draw_CH_chain_simdata <- function(n_draws, sigma_divide = 3, n_cores = 1){

  sims <- sims_small # Change later to use the full set of municipalities

  # Use simulated dataset
  land_parameters_list <- sims$land_params
  V_total_obs <- sims$V_total_obs
  V_cum_obs <- sims$V_cum_obs
  D_flow_obs <- sims$D_flow_obs

  # List of those parameters passed to moment_condition that will not change
  # between iterations of the MCMC run
  fixed_params_list <- list(n_cores = n_cores,
                            land_parameters_list = land_parameters_list,
                            V_total_obs = V_total_obs,
                            V_cum_obs = V_cum_obs,
                            D_flow_obs = D_flow_obs)

  # Starting values, proposal variances, and bounds for uniform priors
  true_params <- unlist(sims$model_params)
  lower <- true_params / c(10, 1, 1, 1, 10, 1, 10, 10)
  upper <- true_params * c(10, 1, 1, 1, 10, 1, 10, 10)
  starting_values <- true_params + c(0.8, 0, 0, 0, 2, 0, 0.75, 0.005)
  sigma <- (upper - lower) / sigma_divide
  n_params <- length(sigma)

  # Matrix of NAs to store the MCMC draws, vector of NAs for the pseudo-Likelihood
  draws <- matrix(NA_real_, nrow = n_draws + 1, ncol = n_params)
  colnames(draws) <- names(sims$model_params)
  Like <- rep(NA, n_draws)

  # Evaluate pseudo-Likelihood at MCMC starting values
  draws[1,] <- starting_values
  MC_temp <- do.call(moment_condition,
                     c(as.list(starting_values), fixed_params_list))
  Like[1] <- Like_CUE(MC_temp)

  # Random walk Metroplis-Hastings updates with normal proposals
  accept_count <- 0L
  for(i in 2:(n_draws + 1)){
    proposal <- draws[i-1,] + rnorm(n_params, mean = 0, sd = sigma)

    # Reject a proposal that lies outside the support of the prior
    if(any(proposal > upper) || any(proposal < lower)){
      draws[i,] <- draws[i-1,]
      Like[i] <- Like[i-1]

    # Otherwise evaluate the pseudo-Likelihood of the proposal
    } else {
      MC_temp <- do.call(moment_condition,
                     c(as.list(proposal), fixed_params_list))
      Like_proposal <- Like_CUE(MC_temp)

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
               upper = upper)
  out <- list(draws = draws, Like = Like, accept_rate = accept_count / n_draws,
              call = call)
  return(out)
}



