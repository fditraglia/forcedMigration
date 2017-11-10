library(forcedMigration)

set.seed(1234)

# Set parameters of simulation
model <- list(delta = 0.2,
              tau_ell = 0.45,
              tau_n = 0.45,
              r = 0.1,
              a0 = 1,
              a1 = 1,
              rho = 0.3)
obs <- list(log_dbar = log(0.001),
            nu = seq(-0.3, 0.3, length.out = 5),
            eta = seq(-0.1, 0.1, length.out = 16)) # One fewer than T
sim_params <- list(model = model, obs = obs)
rm(model, obs)


# Solve model at given values for structural parameters
get_migration_flow_i <- function(i) {
  Vcum_i <- as.matrix(Vcum[i,])
  land_i <- land_parameters[i,]
  par <- sim_params$model
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


# test run
get_migration_flow_i(275)

dstar_list <- lapply(1:nrow(Vcum), get_migration_flow_i)
dstar <- do.call(rbind, dstar_list)
colnames(dstar) <- colnames(Vcum)
rownames(dstar) <- rownames(Vcum)
dstar_lag <- cbind(rep(0, nrow(dstar)), dstar[,-ncol(dstar)])
dstar_friction <- with(sim_params$model, (1 - rho) * dstar + rho * dstar_lag)

# Clean up: only keep fam, pop, and dstar_friction
rm(dstar, dstar_list, dstar_lag, get_migration_flow_i)


# Simulate from observation model
pop <- cross_section$popn1993
mu <- with(sim_params$obs,
           (pop * t(t(dstar_friction + exp(log_dbar)) *
                              exp(c(0, eta)))) %o% exp(nu))
simZ <- array(rpois(length(mu), lambda = mu), dim = dim(mu))

devtools::use_data(simZ, overwrite = TRUE)
devtools::use_data(sim_params, overwrite = TRUE)
rm(list = ls())
