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
w <- c(0.6, 0.3, 0.1)
a <- c(0, 0.2, 0.7)
b <- 0.8
get_tanh_fixed_point(w, a, b)
get_tanh_fixed_point_cpp(w, a, b)
*/

