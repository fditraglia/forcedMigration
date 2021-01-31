#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericVector get_tanh_fixed_point_cpp(NumericVector w, NumericVector a, double b,
                                       double x) {
  return w * tanh(a + b * x);
}

//  x_prev <- start
//  x <- sum(w * tanh(a + b * x_prev))
//  iter_count <- 1L
//
//  while((iter_count <= max_iter) && (abs(x - x_prev) > tol)) {
//    x_prev <- x
//    x <- sum(w * tanh(a + b * x_prev))
//    iter_count <- iter_count + 1L
//  } # END while
//  return(x)

//double get_migration_eq(double V, double start, double delta, double tau_ell,
//                        double tau_n, double r, double a0,
//                        double a1, double p, double q, double H_bar,
//                        double omega_n)
//{
//  const int max_iter = 500;
//  const double tol = 0.0001;
//  int i = 0;
//  double x = start; // Where to start the iteration? Use 0.0 unless you have a bound
//  double f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
//                       omega_n);
//  while((i < max_iter) & (std::abs(f - x) > tol)){
//    x = f;
//    f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
//    i++;
//  }
//  return(f);
//}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
w <- c(0.7, 0.2, 0.1)
a <- c(0, 0.1, 0.5)
b <- 1
get_tanh_fixed_point_cpp(w, a, b, 0.3)
w * tanh(a + b * 0.3)
*/

