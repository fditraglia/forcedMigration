#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
List get_S(arma::mat X_star, arma::vec y_star, arma::mat Q, arma::vec gamma) {
  int n_gamma = gamma.n_elem;

  // Null model
  arma::vec coef0 = arma::solve(X_star, y_star);
  arma::vec resid0 = y_star - X_star * coef0;
  double S0 = arma::dot(resid0, resid0);

  // Alternative model for sequence of potential gammas
  arma::mat coef1(3, n_gamma, arma::fill::zeros);
  arma::vec S1(n_gamma, arma::fill::zeros);

  for(int i = 0; i < n_gamma; ++i){
    arma::mat Z = arma::conv_to<arma::mat>::from(Q > gamma(i));
    arma::vec z_star = arma::vectorise(Z.each_col() - arma::mean(Z, 1));
    arma::mat design_matrix = arma::join_horiz(z_star, X_star);
    coef1.col(i) = arma::solve(design_matrix, y_star);
    arma::vec resid1 = y_star - design_matrix * coef1.col(i);
    S1(i) = arma::dot(resid1, resid1);
  }

  arma::uword gamma_hat_index = S1.index_min();
  double gamma_hat = gamma(gamma_hat_index);
  double S1_gamma_hat = S1(gamma_hat_index);

  arma::vec coef1_gamma_hat = coef1.col(gamma_hat_index);
  arma::mat Z_gamma_hat = arma::conv_to<arma::mat>::from(Q > gamma_hat);
  arma::vec z_star_gamma_hat = arma::vectorise(Z_gamma_hat.each_col() -
    arma::mean(Z_gamma_hat, 1));
  arma::mat design_matrix_gamma_hat = arma::join_horiz(z_star_gamma_hat, X_star);
  arma::vec e_gamma_hat = y_star - design_matrix_gamma_hat * coef1_gamma_hat;

  return List::create(Named("gamma_hat") = gamma_hat,
                      Named("S0") = S0, Named("coef0") = coef0,
                      Named("S1_gamma_hat") = S1_gamma_hat,
                      Named("S1") = S1,
                      Named("coef1_gamma_hat") = coef1_gamma_hat,
                      Named("e_tilde") = resid0,
                      Named("e_gamma_hat") = e_gamma_hat);
}
