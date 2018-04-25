get_migration_deriv_i <- function(i, est, mymethod = 'simple') {
  f <- function(params){
    parlist <- list(delta = params[1],
                    tau_ell = params[2],
                    tau_n = params[3],
                    r = params[4],
                    a0 = params[5],
                    a1 = params[6])
    get_migration_flow_i(i, parlist)
  }
  numDeriv::jacobian(f, est, method = mymethod)
}

get_avar_migration <- function(theta_est, X, ncores = 1) {
  # theta_est is the *full* parameter value, assumes that covariates are included
  # and that we have estimated a model *with* time effects.
  
  fitted_model <- negloglike_outer_D(theta_est, obsZ, X, time_effects = TRUE,
                                     return_inner = TRUE, ncores)
  D_dstar <- parallel::mclapply(1:nrow(cross_section), function(i)
    get_migration_deriv_i(i, as.numeric(fitted_model$parvecs[i,])), mc.cores = ncores)
  
  jBias <- exp(fitted_model$nu)
  tBias <- exp(c(0, fitted_model$eta))
  mu <- t(t(fitted_model$D) * tBias) %o% jBias
  pop <- cross_section$popn1993
  D_mu_dbar <- (pop %o% tBias) %o% jBias
  D_mu_rho <- ((pop %o% tBias) * with(fitted_model, dstar_lag - dstar)) %o% jBias
  
  
  rho <- fitted_model$rho
  nI <- dim(mu)[1]
  nT <- dim(mu)[2]
  nJ <- dim(mu)[3]
  nTheta <- ncol(D_dstar[[1]])
  D_dstar_diff <- lapply(1:length(D_dstar), function(i) (1 - rho) * D_dstar[[i]] +
                           rho * rbind(rep(0, nTheta), D_dstar[[i]][-nT,]))
  D_dstar_diff <- array(as.numeric(unlist(D_dstar_diff)), dim = c(nT, nTheta, nI))
  D_dstar_diff <- aperm(D_dstar_diff, c(3, 1, 2))
  D_dstar_diff <- lapply(1:dim(D_dstar_diff)[3], function(index) D_dstar_diff[,,index])
  names(D_dstar_diff) <- names(X)
  
  f <- function(index) {
    out <- as.matrix(cbind(rep(1, nrow(covariates)), covariates[,X[[index]]]))
    colnames(out) <- c(paste0('const_', names(X)[index]), X[[index]])
    return(out)
  }
  Xvalues <- lapply(1:length(X), f)
  names(Xvalues) <- names(X)
  
  total <- 0
  for(i in 1:nI){
    for(tee in 1:nT) {
      for(j in 1:nJ){
        # Since a1 and tau_ell are *fixed* we don't generate standard errors for these
        D_beta <- D_mu_dbar[i,tee,j] * 
          c(D_dstar_diff$delta[i,tee] * Xvalues$delta[i,] * fitted_model$parvecs[i,]$delta, 
            D_dstar_diff$tau_n[i,tee] * Xvalues$tau_n[i,]* fitted_model$parvecs[i,]$tau_n, 
            D_dstar_diff$r[i,tee] * Xvalues$r[i,] * fitted_model$parvecs[i,]$r, 
            D_dstar_diff$a0[i,tee] * Xvalues$a0[i,] * fitted_model$parvecs[i,]$a0)
        
        D_dbar <- D_mu_dbar[i,tee,j]
        D_rho <- D_mu_rho[i,tee,j]
        D_nu <- rep(0, length(jBias))
        D_nu[j] <- mu[i,tee,j]
        D_eta <- rep(0, length(tBias) - 1)
        if(tee > 1) D_eta[tee - 1] <- mu[i,tee,j]
        D_everything <- c(D_beta, dbar = D_dbar, rho = D_rho, D_nu, D_eta)
        temp <- D_everything %o% D_everything / mu[i,tee,j]
        total <- total + temp
      }
    }
  }
  
  std_residuals <- (obsZ - mu) / sqrt(mu)
  s_hat_sq <- var(std_residuals)
  AVAR <- s_hat_sq * solve(total)
  return(AVAR)
}

