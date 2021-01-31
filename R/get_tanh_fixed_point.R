get_tanh_fixed_point <- function(w, a, b, start = -1, tol = 1e-5, max_iter = 100) {
# Calculate the smallest fixed point of: f(x) = sum(w * tanh(a + b * x))
  stopifnot(is.atomic(w) && is.atomic(a) && (length(w) == length(a)))
  stopifnot(is.atomic(b) && length(b) == 1L)

  x_prev <- start
  x <- sum(w * tanh(a + b * x_prev))
  iter_count <- 1L

  while((iter_count <= max_iter) && (abs(x - x_prev) > tol)) {
    x_prev <- x
    x <- sum(w * tanh(a + b * x_prev))
    iter_count <- iter_count + 1L
  } # END while
  return(x)
} # END get_tanh_fixed_point



plot_tanh_fixed_point <- function(w, a, b) {
# Plot f(x) = sum(w * tanh(a + b * x)) for x = [-1, 1] alongside the 45-degree line
  stopifnot(is.atomic(w) && is.atomic(a) && (length(w) == length(a)))
  stopifnot(is.atomic(b) && length(b) == 1L)
  x <- seq(-1, 1, by = 0.001)
  bX <- outer(b * x, rep.int(1L, length(a)))
  A <- outer(rep.int(1L, length(x)), a)
  y <- drop(tcrossprod(w, tanh(A + bX)))
  plot(x, y, type = 'l', lwd = 2, col = 'blue', xlab = 'x', ylab = 'f(x)')
  abline(0, 1)
}


#------------- Some testing/debugging code
#library(microbenchmark)
#microbenchmark(get_tanh_fixed_point(w = c(0.7, 0.2, 0.1),
#                     a = c(0, 0.1, 0.5),
#                     b = 1))

#test_me <- list(w = c(0.7, 0.2, 0.1),
#                a = c(0, 0.1, 0.5),
#                b = 1)

#do.call(plot_tanh_fixed_point, test_me)
#fixed_point <- do.call(get_tanh_fixed_point, test_me)
#points(fixed_point, fixed_point, pch = 20, cex = 1.5, col = 'red')
