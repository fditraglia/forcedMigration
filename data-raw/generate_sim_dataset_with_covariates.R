library(forcedMigration)

set.seed(1234)

#----------------------- Pump-priming run
get_migration_flow_i(275, list(delta = 10000, tau_ell = 1, tau_n = 0.3,
                               r = 0.8, a0 = 2, a1 = 0))

#----------------------- Set up parameter heterogeneity

# Covariates for each of the model parameters
X <- list(delta = c('bureaucracy', 'offices'),
          tau_ell = NULL, # Intercept only
          tau_n = NULL, # Intercept only
          r = c('rainfall', 'land_return'),
          a0 = c('radio', 'dist_cap', 'elevation'),
          a1 = NULL) # Intercept only

# Calibrated so delta can vary between 4966 and 20138 across municipalities
b_delta <- c(log(10000), # Intercept: "median" municipality has delta = 10000
             -0.7, # Slope for bureaucracy
             -0.7) # Slope for offices

b_tau_ell <- log(1) # Set tau_ell = 1
b_tau_n <- log(0.3) # Set tau_n = 0.3

# Calibrated so r can vary between 0.6 and 1
b_r <- c(log(0.8), # Intercept: "median" municipality has r = 0.8
         0.3, # Slope for rainfall
         0.3) # Slope for land_return

# Calibrated so a0 can vary between 0.74 and 5.44
b_a0 <- c(log(2), # Intercept: "median" municipality has a0 = 2
          2/3, # Slope for radio
          2/3, # Slope for dist_cap
          2/3) # Slope for elevation

b_a1 <- -Inf # Set a1 = 0

# Full set of betas
b_displacement <- c(b_delta, b_tau_ell, b_tau_n, b_r, b_a0, b_a1)

#------------------------ Solve structural model
npars <- lapply(X, function(x) length(x) + 1)
stopifnot(sum(unlist(npars)) == length(b_displacement))

# Indices for elements of b_db_displacement that correspond to a given element of X
par_indices <- lapply(X, function(x) 1:(length(x) + 1))
for(j in 2:length(par_indices)) {
  par_indices[[j]] <- par_indices[[j]] + max(par_indices[[j - 1]])
}
rm(j, npars)
ones <- rep(1, nrow(covariates))
get_par_j <- function(j){
  covariates_j <- as.matrix(cbind(ones, covariates[,X[[j]]]))
  coefficients_j <- b_displacement[par_indices[[j]]]
  return(exp(covariates_j %*% coefficients_j))
}

par_D_outer <- lapply(1:length(X), get_par_j)
par_D_outer <- as.data.frame(par_D_outer)
names(par_D_outer) <- names(X)

solve_model_i <- function(i) {
  par_model_i <- list(delta = par_D_outer$delta[i],
                      tau_ell = par_D_outer$tau_ell[i],
                      tau_n = par_D_outer$tau_n[i],
                      r = par_D_outer$r[i],
                      a0 = par_D_outer$a0[i],
                      a1 = par_D_outer$a1[i])
  get_migration_flow_i(i, par_model_i)
}

dstar <- do.call(rbind, parallel::mclapply(1:nrow(Vcum_pop),
                                             solve_model_i, mc.cores = 8))
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

##----------------------- Solve for payoffs to contract (Dstar, Xstar, Dmax, Xmax)
#payoffs <- do.call(get_payoffs_list, c(par_D_outer, n_cores = 8))
#
##----------------------- Solve for lambda_star and S_e_star at true contract parameters
#par_V_outer <- list(gamma = 0.2, alpha = 0.5)
#contracts <- do.call(get_contracts, c(par_V_outer, list(payoffs = payoffs)))
#
##----------------------- Simulate fixed costs of contract (Delta)
#covariates <- mvtnorm::rmvnorm(nrow(contracts),
#                               sigma = toeplitz(c(1, seq(-0.15, 0.15, 0.05))))
#
## Note: our convention is that Delta = S_e_star - x'beta
#X <- cbind(contracts$S_e, rep(-1, nrow(contracts)), -1 * covariates)
#scaling <- IQR(contracts$S_e) * sqrt(3) / pi
#beta <- (1 / scaling) *  c(1, -median(contracts$S_e), seq(-1, 1, length.out = ncol(covariates)))
#y <- X %*% beta
#
##------------------------ Generate simulated violence from optimal contracts
## Function to draw from zero-truncated poisson
#rtpois <- function(N, lambda) {
#  qpois(runif(N, dpois(0, lambda), 1), lambda)
#}
#
#Vstar <- rtpois(nrow(contracts), contracts$lam)
#logit_errors <- rlogis(nrow(contracts))
#contract_indicator <- y > logit_errors
#Vobs <- drop(contract_indicator * Vstar)

#------------------------- Store simulated data and parameters
simulation_covariates <- list(par_D_inner = par_D_inner,
                              par_D_outer = par_D_outer,
                              #par_V_outer = par_V_outer,
                              #beta = beta,
                              b_displacement = b_displacement,
                              dstar = dstar,
                              dstar_lag = dstar_lag,
                              true_displacement = true_displacement,
                              covariates = X,
                              Z = Z)
                   #payoffs = payoffs,
                   #contracts = contracts,
                   #covariates = covariates,
                   #X = X,
                   #y = y,
                   #Vstar = Vstar,
                   #logit_errors = logit_errors,
                   #contract_indicator = contract_indicator,
                   #Vobs = Vobs


devtools::use_data(simulation_covariates, overwrite = TRUE)

#----------------------- Clean up
rm(list = ls())
