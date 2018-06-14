# Function to construct heterogenous parameters of the form exp(x'\theta)
get_heterog_pars <- function(theta, X) {
  # theta is a vector of "common parameters" that are used to construct the
  # heterogeneous parameters via heterog_pars = exp(x ' theta)

  # X is a list that specifies which elements of the data.frame `covariates`
  # are used to construct the heterogenous parameters. For example, on the
  # migration side of the problem, we might have
  #
  #         X <- list(delta = c('bureaucracy', 'offices'),
  #                   tau_ell = NULL,
  #                   tau_n = NULL,
  #                   r = c('rainfall', 'land_return'),
  #                   a0 = c('radio', 'dist_cap', 'elevation'),
  #                   a1 = NULL)
  #
  # where an entry of NULL indicates a parameter that will not vary with x
  # and hence will only have an intercept term associated with it in x' \theta

  # get_heterog_pars parses the list X, creates the appropriate "design matrix"
  # with columns of ones to serve as intercepts, and then calculates exp(x'\theta)

  # Ensure that the length of the parameter agrees with X (adding intercepts)
  npars <- lapply(X, function(x) length(x) + 1)
  stopifnot(sum(unlist(npars)) == length(theta))

  # Indices for elements of theta that correspond to a given element of X
  par_indices <- lapply(X, function(x) 1:(length(x) + 1))
  for(j in 2:length(par_indices)) {
    par_indices[[j]] <- par_indices[[j]] + max(par_indices[[j - 1]])
  }
  rm(j, npars)

  ones <- rep(1, nrow(covariates))
  get_par_j <- function(j){
    covariates_j <- as.matrix(cbind(ones, covariates[,X[[j]]]))
    coefficients_j <- theta[par_indices[[j]]]
    return(exp(covariates_j %*% coefficients_j))
  }
  heterog_pars <- lapply(1:length(X), get_par_j)
  heterog_pars <- as.data.frame(heterog_pars)
  names(heterog_pars) <- names(X)
  return(heterog_pars)
}

