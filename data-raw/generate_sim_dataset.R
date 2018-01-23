library(forcedMigration)

set.seed(1234)

#----------------------- Pump-priming run
get_migration_flow_i(275, list(delta = 10000, tau_ell = 1, tau_n = 0.3,
                               r = 0.8, a0 = 2, a1 = 0))

#----------------------- Solve for eq. migration given observed violence flows
par_D_outer <- list(delta = 10000,
                tau_ell = 1,
                tau_n = 0.3,
                r = 0.8,
                a0 = 2,
                a1 = 0)

dstar <- do.call(rbind, lapply(1:nrow(Vcum),
                               function(i) get_migration_flow_i(i, par_D_outer)))
colnames(dstar) <- colnames(Vcum)
rownames(dstar) <- rownames(Vcum)

#----------------------- Simulate observed displacement measures
par_D_inner <- list(dbar = 0.001,
                rho = 0.3,
                nu = seq(-0.3, 0.3, length.out = 5),
                eta = seq(-0.1, 0.1, length.out = 16)) # One fewer than T

dstar_lag <- cbind(rep(0, nrow(dstar)), dstar[,-ncol(dstar)])
true_displacement <- with(par_D_inner, dbar + (1 - rho) * dstar + rho * dstar_lag)

mu <- with(par_D_inner, (cross_section$popn1993 * t(t(true_displacement) *
                                                      exp(c(0, eta)))) %o% exp(nu))
Z <- array(rpois(length(mu), lambda = mu), dim = dim(mu))

#----------------------- Solve for payoffs to contract (Dstar, Xstar, Dmax, Xmax)
payoffs <- do.call(get_payoffs_list, c(par_D_outer, n_cores = 8))

#----------------------- Solve for lambda_star and S_e_star at true contract parameters
par_V_outer <- list(gamma = 0.2, alpha = 0.5)
contracts <- do.call(get_contracts, c(par_V_outer, list(payoffs = payoffs)))

#----------------------- Simulate fixed costs of contract (Delta)
covariates <- mvtnorm::rmvnorm(nrow(contracts),
                               sigma = toeplitz(c(1, seq(-0.15, 0.15, 0.05))))

# Note: our convention is that Delta = S_e_star - x'beta
X <- cbind(contracts$S_e, rep(-1, nrow(contracts)), -1 * covariates)
scaling <- IQR(contracts$S_e) * sqrt(3) / pi
beta <- (1 / scaling) *  c(1, -median(contracts$S_e), seq(-1, 1, length.out = ncol(covariates)))
y <- X %*% beta

#------------------------ Generate simulated violence from optimal contracts
# Function to draw from zero-truncated poisson
rtpois <- function(N, lambda) {
  qpois(runif(N, dpois(0, lambda), 1), lambda)
}

Vstar <- rtpois(nrow(contracts), contracts$lam)
logit_errors <- rlogis(nrow(contracts))
contract_indicator <- y > logit_errors
Vobs <- drop(contract_indicator * Vstar)

#------------------------- Store simulated data and parameters
simulation <- list(par_D_inner = par_D_inner,
                   par_D_outer = par_D_outer,
                   par_V_outer = par_V_outer,
                   beta = beta,
                   dstar = dstar,
                   dstar_lag = dstar_lag,
                   true_displacement = true_displacement,
                   Z = Z,
                   payoffs = payoffs,
                   contracts = contracts,
                   covariates = covariates,
                   X = X,
                   y = y,
                   Vstar = Vstar,
                   logit_errors = logit_errors,
                   contract_indicator = contract_indicator,
                   Vobs = Vobs)

devtools::use_data(simulation, overwrite = TRUE)

#----------------------- Clean up
rm(list = ls())
