library(forcedMigration)

# set.seed(57)
#
# # Take a small random sample of municipalities
# n_sample <- 48
# municipalities <- cross_section$municipality
# municipality_sample <- sample(municipalities, n_sample)
# small_cross <- subset(cross_section, municipality %in% municipality_sample)
# land_parameters <- small_cross[,c('p', 'q', 'H_bar', 'omega_n')]
# land_parameters_list <- lapply(seq_len(nrow(land_parameters)),
#                                function(i) as.list(land_parameters[i,]))
#
# rm(municipalities, municipality_sample, small_cross, land_parameters)

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

sims <- simulate_from_model(delta, tau_ell, tau_n, r, a0, a1,
                            land_parameters_list, gamma, alpha, n_cores = 8)

sims$land_params <- land_parameters_list
sims$model_params <- list(delta = delta,
                          tau_ell = tau_ell,
                          tau_n = tau_n,
                          r = r,
                          a0 = a0,
                          a1 = a1,
                          gamma = gamma,
                          alpha = alpha)

devtools::use_data(sims, overwrite = TRUE)
rm(list = ls())
