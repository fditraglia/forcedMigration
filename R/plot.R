# params is a list with names that match the arguments of get_migration_eq
plot_migration_eq <- function(params){
  eq <- do.call(get_migration_eq, params)
  x <- seq(0, 1, 0.001)
  fx <- rep(NA_real_, length((x)))
  for(i in 1:length(x)){
    fx[i] <- do.call(get_Dstar, c(D_e = x[i], params))
  }
  plot(x, fx, type = 'l', col = 'blue', xlab = 'D_e', ylab = 'D_star')
  abline(0, 1, col = 'red')
  abline(v = eq, lty = 2, col = 'black')
}

# params is a list with names that match the arguments of get_surplus
plot_surplus <- function(params, V_lower = 0, V_upper = 10, type = 'l'){
  V_seq <- seq(V_lower, V_upper, length.out = 500)
  surplus <- rep(NA_real_, length((V_seq)))
  for(i in 1:length(V_seq)){
    surplus[i] <- do.call(get_surplus, c(V = V_seq[i], params))
  }
  plot(V_seq, surplus, type = type, col = 'blue', xlab = 'V', ylab = 'S(V)')
}

