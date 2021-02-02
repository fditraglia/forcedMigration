#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
double get_tanh_fixed_point_cpp(NumericVector w, NumericVector a, double b,
                                       double start = -1) {
  const double tol = 0.00001;
  const int max_iter = 500;
  int i = 0;
  double x = start;
  double f = sum(w * tanh(a + b * x));

  while((i < max_iter) & (std::abs(f - x) > tol)){
    x = f;
    f = sum(w * tanh(a + b * x));
    i++;
  }
  return f;
}

/*** R
# This example is based on Municipality 13670 (with rounding)
Vcum <- 121
nFam <- 3170
w <- c(0.83, 0.029, 0.0202, 0.0177, 0.0271, 0.0183, 0.0139, 0.0268, 0.0095,
       0.0044, 0.0025, 6e-04, 0, 0)

h <- c(0, 0.35, 1.93, 3.91, 7.62, 12.49, 17.74, 30.25, 71.57, 147.95, 294.79,
       795.24, 0, 0)
tau <- 2.5
r <- 0.05
J <- 1
a <- tau * (Vcum / nFam) - r * h
(1 + sum(w * tanh(a))) / 2 # Fraction who leave with no social interactions (J=0)

plot_tanh_fixed_point(w, a, J)
get_tanh_fixed_point(w, a, J)
get_tanh_fixed_point_cpp(w, a, J)
microbenchmark::microbenchmark(get_tanh_fixed_point(w, a, J),
                               get_tanh_fixed_point_cpp(w, a, J))
*/

