// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <Rcpp.h>
#include <RcppNumerical.h>
#include "brent.hpp"

using namespace Numer;
using namespace Rcpp;

double get_Q(double V, double D_e, double delta){
  const double V_tilde = V / (1 - D_e);
  if(fabs(1.0 - delta) > 0.00001){
    return((delta / (1 - delta)) * (1 / (delta + (1 - delta) * exp(-V_tilde)) - 1));
  }else{
    return(1 - exp(-V_tilde));
  }
}

class MigrationIntegrand: public Func
{
  private:
    double Q, tau_ell, r, a0, a1, p, q, H_bar;
  public:
    MigrationIntegrand(double Q_, double tau_ell_,
                       double r_, double a0_,
                       double a1_, double p_,
                       double q_, double H_bar_): Q(Q_), tau_ell(tau_ell_),
                                                  r(r_), a0(a0_), a1(a1_),
                                                  p(p_), q(q_), H_bar(H_bar_) {}

    double operator()(const double& h) const
    {
      double mu = H_bar * p / (p + q);
      double upper = (exp(tau_ell * Q - r * h / mu) - 1);
      double log_F_c_given_h = R::pbeta(upper, a0 + a1 * h / mu, 1, 1, 1);
      double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
      return exp(log_F_c_given_h + log_f_h);
    }
};

// [[Rcpp::export]]
double get_Dstar(double D_e, double V, double delta, double tau_ell,
                 double tau_n, double r, double a0, double a1,
                 double p, double q, double H_bar, double omega_n)
{
  if(V > 0.0){
    // Dis-utility of violence
    const double Q = get_Q(V, D_e, delta);
    const double mu = H_bar * p / (p + q);

    // Landholding families
    MigrationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
    double err_est;
    int err_code;
    const double D_ell = integrate(g, 0.0, mu * tau_ell * Q / r, err_est, err_code);

    // Landless families
    const double D_n = R::pbeta(exp(tau_n * Q) - 1, a0, 1, 1, 0);

    // Overall migration
    return omega_n * D_n + (1 - omega_n) * D_ell;

  } else {
    return 0.0;
  }
}


// [[Rcpp::export]]
double get_migration_eq(double V, double start, double delta, double tau_ell,
                        double tau_n, double r, double a0,
                        double a1, double p, double q, double H_bar,
                        double omega_n)
{
  const int max_iter = 500;
  const double tol = 0.0001;
  int i = 0;
  double x = start; // Where to start the iteration? Use 0.0 unless you have a bound
  double f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                       omega_n);
  while((i < max_iter) & (std::abs(f - x) > tol)){
    x = f;
    f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
    i++;
  }
  return(x);
}

// [[Rcpp::export]]
NumericVector get_migration_cum(NumericVector V_cum, double delta, double tau_ell,
                                double tau_n, double r, double a0, double a1,
                                double p, double q, double H_bar, double omega_n)
{
  int n = V_cum.length();
  NumericVector Dstar_cum(n);
  double start = 0.0;
  double Dstar_temp;
  for(int i = 0; i < n; ++i) {
    Dstar_temp = get_migration_eq(V_cum[i], start, delta, tau_ell, tau_n, r,
                                  a0, a1, p, q, H_bar, omega_n);
    Dstar_cum[i] = Dstar_temp;
    start = Dstar_temp;
  }
  return Dstar_cum;
}


class ExpropriationIntegrand: public Func
{
  private:
    double Q, tau_ell, r, a0, a1, p, q, H_bar;
  public:
    ExpropriationIntegrand(double Q_, double tau_ell_, double r_,
                           double a0_, double a1_, double p_, double q_,
                           double H_bar_): Q(Q_), tau_ell(tau_ell_),
                                           r(r_), a0(a0_), a1(a1_),
                                           p(p_), q(q_), H_bar(H_bar_) {}
    double operator()(const double& h) const
    {
      const double mu = H_bar * p / (p + q);
      double upper = (exp(tau_ell * Q - r * h / mu) - 1);
      double log_F_c_given_h = R::pbeta(upper, a0 + a1 * h / mu, 1, 1, 1);
      double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
      return exp(log(h) + log_F_c_given_h + log_f_h);
    }
};

// [[Rcpp::export]]
double get_surplus(double V, double delta, double tau_ell, double tau_n,
                   double r, double a0, double a1, double p, double q,
                   double H_bar, double omega_n, double gamma, double alpha)
{
  if(V > 0.0){
    // Equilibrium migration at violence level V
    double Dstar = get_migration_eq(V, 0.0, delta, tau_ell, tau_n, r, a0,
                          a1, p, q, H_bar, omega_n);
    // Equilibrium dis-utility of violence at violence level V
    const double Q = get_Q(V, Dstar, delta);
    // Calculate total land expropriated (in per capita terms)
    const double mu = H_bar * p / (p + q);
    ExpropriationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
    double err_est;
    int err_code;
    const double res = integrate(g, 0.0, mu * tau_ell * Q / r, err_est, err_code);
    double Xstar = (1 - omega_n) * res;
    return Xstar - gamma * Dstar - alpha * V / (1 - Dstar);
  } else {
    return 0.0;
  }
}

// [[Rcpp::export]]
double get_surplus_infeas(double V_tilde, double delta, double tau_ell,
                          double tau_n, double r, double a0, double a1,
                          double p, double q, double H_bar, double omega_n,
                          double gamma, double alpha)
{
// This is an "infeasible" surplus function that pretends we could freely
// choose V_tilde = V / (1 - Dstar). The idea is as follows. Imagine that a
// dictator could freely assign Q, the dis-utility of violence. Q is a function
// of V_tilde, but in equilibrium some values of V_tilde cannot be obtained if
// there is tipping. When this occurs, certain values of Q cannot be attained
// and hence certain levels of migration cannot be attained. Nevertheless,
// because households decisions depend only on Q, it remains a well-defined
// question to ask what they *would* have at any particular value of Q.
// To carry out this exercise, we will not solve a fixed point problem
// to determine equilibrium migration. Instead we will set Q via V_tilde. An
// easy way to re-use our code from above rather than writing a new function
// Q(V_tilde) is simply to call Q(V_tilde, D_e = 0, delta). Once we have the
// value of Q, we need to know who would migrate, and how much land they would
// leave behind. This involves integrating both the expropriation integrand
// and the migration integrand. But again: it does not involve imposing
// equilibrium and hence there is no fixed point problem to solve.

  if(V_tilde > 0.0){
    // Potentially infeasible Q that results from V_tilde, which we obtain by
    // setting D_e = 0.0 and V = V_tilde
    const double Q = get_Q(V_tilde, 0.0, delta);

    // Construct expropriation integrand using the potentially infeasible
    // level of migration D_infeas
    ExpropriationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
    double err_est;
    int err_code;
    const double mu = H_bar * p / (p + q);
    const double res = integrate(g, 0.0, mu * tau_ell * Q / r, err_est, err_code);
    double Xstar = (1 - omega_n) * res;

    // Potentially infeasible migration from the Q that results from V_tilde,
    // which we obtain by setting D_e = 0.0 and V = V_tilde
    const double D_infeas = get_Dstar(0.0, V_tilde, delta, tau_ell, tau_n, r,
                                      a0, a1, p, q, H_bar, omega_n);
    return Xstar - gamma * D_infeas - alpha * V_tilde;
  } else {
    return 0.0;
  }
}

// [[Rcpp::export]]
double get_X_max(double tau_ell, double r, double a0, double a1,
                 double p, double q, double H_bar, double omega_n)
{
// This function calculates the maximum amount of land (per capita) that can
// be expropriated, given model parameters. This is essentially a limit as
// violence V (or infeasible violence V_tilde ) goes to infinity. We calculate
// this by setting Q = 1.0 when setting up the expropriation integrand.
  const double Q = 1.0;
  ExpropriationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  const double mu = H_bar * p / (p + q);
  const double res = integrate(g, 0.0, mu * tau_ell * Q / r, err_est, err_code);
  return (1 - omega_n) * res;
}


// [[Rcpp::export]]
double get_D_max(double tau_ell, double tau_n, double r, double a0, double a1,
                 double p, double q, double H_bar, double omega_n)
{
// This function calculates the maximum fraction of families that can be
// displaced, given model parameters. This is essentially a limit as
// violence V (or infeasible violence V_tilde ) goes to infinity. We calculate
// this by setting Q = 1.0 when setting up the migration integrand.
  const double Q = 1.0;
  const double mu = H_bar * p / (p + q);

  // Landholding families
  MigrationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  const double D_ell = integrate(g, 0.0, mu * tau_ell * Q / r, err_est, err_code);

  // Landless families
  const double D_n = R::pbeta(exp(tau_n * Q) - 1, a0, 1, 1, 0);

  // Overall migration
  return omega_n * D_n + (1 - omega_n) * D_ell;
}

// Functor to pass to the BRENT optimization routine to maximize infeasible
// surplus. Note that BRENT *minimizes* so the functor calculates *negative*
// surplus.
class negInfeasibleSurplus : public brent::func_base
{
private:
  double delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n, gamma, alpha;
public:
  negInfeasibleSurplus (double delta_, double tau_ell_, double tau_n_,
                        double r_, double a0_, double a1_, double p_,
                        double q_, double H_bar_, double omega_n_,
                        double gamma_, double alpha_) : delta(delta_),
                                                        tau_ell(tau_ell_),
                                                        tau_n(tau_n_),
                                                        r(r_), a0(a0_), a1(a1_),
                                                        p(p_), q(q_),
                                                        H_bar(H_bar_),
                                                        omega_n(omega_n_),
                                                        gamma(gamma_),
                                                        alpha(alpha_) {}
  double operator() (double V_tilde) {
    return -1.0 * get_surplus_infeas(V_tilde, delta, tau_ell, tau_n, r, a0, a1,
                                     p, q, H_bar, omega_n, gamma, alpha);
  }
};

// [[Rcpp::export]]
double get_V_tilde_star(double delta, double tau_ell, double tau_n,
                        double r, double a0, double a1, double p, double q,
                        double H_bar, double omega_n, double gamma,
                        double alpha){

  negInfeasibleSurplus neg_V_tilde(delta, tau_ell, tau_n, r, a0, a1, p, q,
                                   H_bar, omega_n, gamma, alpha);
  double X_max = get_X_max(tau_ell, r, a0, a1, p, q, H_bar, omega_n);
  double upper = X_max / alpha;
  double lower = 0.0;
  double V_tilde_star;
  double neg_S_tilde_star = brent::local_min(lower, upper, 0.0001,
                                             neg_V_tilde, V_tilde_star);
  return V_tilde_star;
}


// Object to contain result of bisection algorithm
struct bisectObj{
  double x_lower, x_upper;
  double y_lower, y_upper;
  int n_iter;
};

bisectObj get_V(double V_tilde, double delta, double tau_ell, double tau_n,
               double r, double a0, double a1, double p, double q, double H_bar,
               double omega_n){

  double tol = 0.005;
  double D_max = get_D_max(tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
  double V_L = V_tilde * (1.0 - D_max) - 0.1 * tol;
  double V_U = V_tilde + 0.1 * tol;
  int n_test = 10;
  double inc = (V_U - V_L) / double(n_test - 1);
  IntegerVector seq = seq_len(n_test) - 1;
  NumericVector V_test = V_L + inc * as<NumericVector>(seq);
  NumericVector D_test = get_migration_cum(V_test, delta, tau_ell, tau_n,
                                           r, a0, a1, p, q, H_bar, omega_n);
  NumericVector V_tilde_test = V_test / (1 - D_test);

  //Which element of V_tilde_test is the *first* that exceeds V_tilde?
  int i = 0;
  while(V_tilde_test[i] < V_tilde) ++i;
  V_L = V_test[i - 1];
  V_U = V_test[i];
  double start = D_test[i - 1];
  double V_tilde_L = V_L;
  double V_tilde_U = V_U;
  int n_iter = 0;

  while(std::abs(V_U - V_L) > tol && (n_iter <= 50)){
    double V_M = 0.5 * (V_L + V_U);
    double Dstar_M = get_migration_eq(V_M, start, delta, tau_ell, tau_n, r,
                                      a0, a1, p, q, H_bar, omega_n);
    double V_tilde_M = V_M / (1 - Dstar_M);
      if(V_tilde_M < V_tilde){
        V_L = V_M;
        V_tilde_L = V_tilde_M;
        start = Dstar_M;
      } else {
        V_U = V_M;
        V_tilde_U = V_tilde_M;
      }
      n_iter += 1;
  }
  bisectObj out;
  out.n_iter = n_iter;
  out.x_lower = V_L;
  out.x_upper = V_U;
  out.y_lower = V_tilde_L;
  out.y_upper = V_tilde_U;
  return out;
}

// [[Rcpp::export]]
List get_V_cpp(double V_tilde, double delta, double tau_ell, double tau_n,
               double r, double a0, double a1, double p, double q, double H_bar,
               double omega_n){

  bisectObj out = get_V(V_tilde, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                        omega_n);

  return List::create(Named("V_L") = out.x_lower,
                      Named("V_U") = out.x_upper,
                      Named("V_tilde_L") = out.y_lower,
                      Named("V_tilde_U") = out.y_upper,
                      Named("n_iter") = out.n_iter);
}

// [[Rcpp::export]]
double get_V_star_cpp(double delta, double tau_ell, double tau_n, double r,
                      double a0, double a1, double p, double q, double H_bar,
                      double omega_n, double gamma, double alpha){

  double V_tilde_star = get_V_tilde_star(delta, tau_ell, tau_n, r, a0, a1, p, q,
                                         H_bar, omega_n, gamma, alpha);

  double V_star = 0.0; // The default return value is zero
  double tol = 0.001;
  if(fabs(V_tilde_star) > tol){
    bisectObj V_bisect = get_V(V_tilde_star, delta, tau_ell, tau_n, r, a0, a1,
                             p, q, H_bar, omega_n);
    double S_L = get_surplus_infeas(V_bisect.y_lower, delta, tau_ell, tau_n, r,
                                                      a0, a1, p, q, H_bar,
                                                      omega_n, gamma, alpha);
    double S_U = get_surplus_infeas(V_bisect.y_upper, delta, tau_ell, tau_n, r,
                                                      a0, a1,  p, q, H_bar,
                                                      omega_n, gamma, alpha);
    if((S_L > 0.0) && (S_U > 0.0)){
      if(S_L > S_U) V_star = V_bisect.x_lower;
      else V_star = V_bisect.x_upper;
    }
  }
  return V_star;
}
