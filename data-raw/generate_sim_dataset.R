library(forcedMigration)

land_parameters <- cross_section[, c('p', 'q', 'H_bar', 'omega_n')]
land_parameters_list <- lapply(1:nrow(land_parameters),
                               function(i) as.list(land_parameters[i,]))
rm(land_parameters)

# "Fixed" Model parameters
delta <- 0.2
tau_ell <- 0.6
tau_n <- 0.6
r <- 0.1
a0 <- 5
a1 <- 0
gamma <- 0.25
alpha <- 0.005

set.seed(2619)
sims <- simulate_from_model(delta, tau_ell, tau_n, r, a0, a1,
                            land_parameters_list, gamma, alpha, n_cores = 4)

sims$land_params <- land_parameters_list
sims$model_params <- list(delta = delta,
                          tau_ell = tau_ell,
                          tau_n = tau_n,
                          r = r,
                          a0 = a0,
                          a1 = a1,
                          gamma = gamma,
                          alpha = alpha)

# Take a random sample to create a "small" sim dataset
random_municipalities <- sample(1:length(sims$V_total), size = 100)
sims_small <- list(V_total = sims$V_total[random_municipalities],
                   V_total_obs = sims$V_total_obs[random_municipalities],
                   V_cum_obs = sims$V_cum_obs[random_municipalities],
                   D_flow = sims$D_flow[random_municipalities],
                   D_flow_obs = sims$D_flow_obs[random_municipalities],
                   land_params = sims$land_params[random_municipalities],
                   model_params = sims$model_params)


devtools::use_data(sims, overwrite = TRUE)
devtools::use_data(sims_small, overwrite = TRUE)
rm(list = ls())
