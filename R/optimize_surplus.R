get_V_R <- function(V_tilde, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n){
  tol <- 0.005
  n <- 10
  D_max <- get_D_max(tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n)
  V_test <- seq(V_tilde * (1 - D_max) - tol/10, V_tilde + tol/10, length.out = n)
  D_test <- get_migration_cum(V_test, delta, tau_ell, tau_n, r, a0, a1, p, q,
                              H_bar, omega_n)
  V_tilde_test <- V_test / (1 - D_test)
  L_test <- max(which(V_tilde_test < V_tilde))
  V_L <- V_test[L_test]
  V_U <- V_test[L_test + 1]
  start <- D_test[L_test]
  n_iter <- 0L
  while((abs(V_U - V_L) > tol) && (n_iter <= 50L)){
    V_M <- 0.5 * (V_L + V_U)
    D_M <- get_migration_eq(V_M, start, delta, tau_ell, tau_n, r, a0, a1, p, q,
                            H_bar, omega_n)
    V_tilde_M <- V_M / (1 - D_M)
    if(V_tilde_M < V_tilde){
     V_L <- V_M
     V_tilde_L <- V_tilde_M
     start <- D_M
    } else {
     V_U <- V_M
     V_tilde_U <- V_tilde_M
    }
    n_iter <- n_iter + 1L
  }
  return(list(V_L = V_L, V_U = V_U,
              V_tilde_L = V_tilde_L, V_tilde_U = V_tilde_U,
              n_iter = n_iter))
}


get_V_star_R <- function(delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n,
                       gamma, alpha){

  X_max <- get_X_max(tau_ell, r, a0, a1, p, q, H_bar, omega_n)

  g <- function(V_tilde){
    get_surplus_infeas(V_tilde, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                       omega_n, gamma, alpha)
  }

  V_star_tilde <- optimize(g, lower = 0, upper = X_max / alpha,
                           maximum = TRUE)$maximum
  if(abs(V_star_tilde) > 0.001){
    bisect_result <- get_V_R(V_star_tilde, delta, tau_ell, tau_n, r, a0, a1, p, q,
                           H_bar, omega_n)
    S_L <- get_surplus_infeas(bisect_result$V_tilde_L, delta, tau_ell, tau_n, r,
                              a0, a1, p, q, H_bar, omega_n, gamma, alpha)
    S_U <- get_surplus_infeas(bisect_result$V_tilde_U, delta, tau_ell, tau_n, r,
                              a0, a1,  p, q, H_bar, omega_n, gamma, alpha)
    if((S_L < 0) && (S_U < 0)){
      V_star <- 0.0
    } else {
      if(S_L > S_U){
        V_star <- bisect_result$V_L
      } else {
        V_star <- bisect_result$V_U
      }
    }
  } else {
    V_star <- 0.0
  }
  return(V_star)
}

