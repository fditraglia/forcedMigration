// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Numer;

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
}


// [[Rcpp::export]]
double get_migration_eq(double V, double delta, double tau_ell,
                        double tau_n, double r, double a0,
                        double a1, double p, double q, double H_bar,
                        double omega_n)
{
  const int max_iter = 500;
  const double tol = 0.0001;
  int i = 0;
  double x = 0.0;
  double f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar,
                       omega_n);
  while((i < max_iter) & (std::abs(f - x) > tol)){
    x = f;
    f = get_Dstar(x, V, delta, tau_ell, tau_n, r, a0, a1, p, q, H_bar, omega_n);
    i++;
  }
  return(x);
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
  // Equilibrium migration at violence level V
  double Dstar = get_migration_eq(V, delta, tau_ell, tau_n, r, a0,
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
