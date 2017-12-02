library(forcedMigration)

set.seed(1234)

# Set parameters of simulation
model <- list(delta = 5000,
              tau_ell = 0.5,
              tau_n = 0.3,
              r = 0.05,
              a0 = 1,
              a1 = 0)
obs <- list(dbar = 0.001,
            rho = 0.3,
            nu = seq(-0.3, 0.3, length.out = 5),
            eta = seq(-0.1, 0.1, length.out = 16)) # One fewer than T
sim_params <- list(model = model, obs = obs)
rm(model, obs)

# Pump-priming run
get_migration_flow_i(275, sim_params$model)

# Solve model at given values for structural parameters
f <- function(i) get_migration_flow_i(i, sim_params$model)
sim_dstar <- do.call(rbind, lapply(1:nrow(Vcum), f))
colnames(sim_dstar) <- colnames(Vcum)
rownames(sim_dstar) <- rownames(Vcum)
dstar_lag <- cbind(rep(0, nrow(sim_dstar)), sim_dstar[,-ncol(sim_dstar)])
dstar_friction <- with(sim_params$obs, (1 - rho) * sim_dstar + rho * dstar_lag)

# Simulate from observation model
pop <- cross_section$popn1993
mu <- with(sim_params$obs,
           (pop * t(t(dstar_friction + dbar) * exp(c(0, eta)))) %o% exp(nu))
simZ <- array(rpois(length(mu), lambda = mu), dim = dim(mu))

devtools::use_data(simZ, overwrite = TRUE)
devtools::use_data(sim_params, overwrite = TRUE)
devtools::use_data(sim_dstar, overwrite = TRUE)
rm(list = ls())
