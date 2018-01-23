#include <vector>
#include <RcppNumerical.h>
#include "brent.hpp"

using namespace Numer;
using namespace Rcpp;

// [[Rcpp::export]]
double get_Q(double V, double D_e, double delta){
  const double V_tilde = V / (1.0 - D_e); // V is assumed to be in per-capita terms
  return(1.0 - exp(-delta * V_tilde));
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
      double med = H_bar * R::qbeta(0.5, p, q, 1, 0);
      double upper = (exp(tau_ell * Q - r * h / med) - 1.0);
      if(upper <= 0.0) {
        return(0.0);
      } else {
        double log_F_c_given_h = R::pbeta(upper, 1.0 + a0 + a1 * h / med, 1, 1, 1);
        double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
        return exp(log_F_c_given_h + log_f_h);
      }
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
    double med = H_bar * R::qbeta(0.5, p, q, 1, 0);
    // Landholding families
    MigrationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
    double err_est;
    int err_code;
    const double D_ell = integrate(g, 0.0, med * tau_ell * Q / r, err_est, err_code);

    // Landless families
    const double D_n = R::pbeta(exp(tau_n * Q) - 1.0, 1.0 + a0, 1, 1, 0);

    // Overall migration
    return omega_n * D_n + (1.0 - omega_n) * D_ell;

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
  return(f);
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
    if(i == 0){
      Dstar_cum[i] = Dstar_temp;
      start = Dstar_temp;
    } else {
      Dstar_cum[i] = fmax(Dstar_temp, Dstar_cum[i - 1]);
      start = Dstar_cum[i];
    }
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
      double med = H_bar * R::qbeta(0.5, p, q, 1, 0);
      double upper = (exp(tau_ell * Q - r * h / med) - 1);
      double log_F_c_given_h = R::pbeta(upper, 1.0 + a0 + a1 * h / med, 1, 1, 1);
      double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
      return exp(log(h) + log_F_c_given_h + log_f_h);
    }
};

// [[Rcpp::export]]
double get_X(double Q, double tau_ell, double r, double a0, double a1,
             double p, double q, double H_bar, double omega_n)
{
// This function calculates the maximum amount of land (per capita)
// expropriated when the dis-utility of violence equals Q, given model
// parameters. Q need not be a value that can be supported in equilibrium. If
// Q is an equilibrium value, the result is equilibrium land expropriated.
  ExpropriationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  double med = H_bar * R::qbeta(0.5, p, q, 1, 0);
  const double res = integrate(g, 0.0, med * tau_ell * Q / r, err_est, err_code);
  return (1.0 - omega_n) * res;
}

// [[Rcpp::export]]
double get_X_max(double tau_ell, double r, double a0, double a1,
                 double p, double q, double H_bar, double omega_n)
{
// This function calculates the maximum amount of land (per capita) that can
// be expropriated, given model parameters. This is essentially a limit as
// violence to infinity. We calculate this by setting Q = 1.0 and calling get_X
  return get_X(1.0, tau_ell, r, a0, a1, p, q, H_bar, omega_n);
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
  double med = H_bar * R::qbeta(0.5, p, q, 1, 0);

  // Landholding families
  MigrationIntegrand g(Q, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  const double D_ell = integrate(g, 0.0, med * tau_ell * Q / r, err_est, err_code);

  // Landless families
  const double D_n = R::pbeta(exp(tau_n * Q) - 1.0, 1.0 + a0, 1, 1, 0);

  // Overall migration
  return omega_n * D_n + (1.0 - omega_n) * D_ell;
}

// [[Rcpp::export]]
List get_payoffs(double delta, double tau_ell, double tau_n, double r,
                 double a0, double a1, double p, double q, double H_bar,
                 double omega_n, double popn)
{
  double D_max = get_D_max(tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
  double X_max = get_X_max(tau_ell, r, a0, a1, p, q, H_bar, omega_n);

  std::vector<double> Dstar;
  std::vector<double> Xstar;

  unsigned int V_realized = 1; // Careful about types here!
  double Dstar_temp = get_migration_eq(double(V_realized) / popn, 0.0, delta,
                                       tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                                       omega_n);
  double Qstar_temp = get_Q(double(V_realized) / popn, Dstar_temp, delta);
  double Xstar_temp = get_X(Qstar_temp, tau_ell, r, a0, a1, p, q, H_bar, omega_n);

  Dstar.push_back(Dstar_temp);
  Xstar.push_back(Xstar_temp);

  // Tolerance of 1 person per 10,000 for Dstar relative to D_max and
  // 0.01 hectare per family for Xstar relative to X_max. Note that the
  // tolerance for Dstar should *not* be set smaller than the tolerance
  // for the fixed point iteration in get_migration_eq. Maximum of 100
  // iterations to ensure that the loop terminates.
  while(((std::abs(Dstar_temp - D_max) > 0.0001) ||
        (std::abs(Xstar_temp - X_max) > 0.001)) && (V_realized <= 100))
  {
    ++V_realized;
    Dstar_temp = get_migration_eq(double(V_realized) / popn, Dstar_temp, delta,
                                  tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
    Dstar.push_back(Dstar_temp);
    Qstar_temp = get_Q(double(V_realized) / popn, Dstar_temp, delta);
    Xstar_temp = get_X(Qstar_temp, tau_ell, r, a0, a1, p, q, H_bar, omega_n);
    Xstar.push_back(Xstar_temp);
  }

  return List::create(Named("D_max") = D_max,
                      Named("X_max") = X_max,
                      Named("Dstar") = Dstar,
                      Named("Xstar") = Xstar);
}
