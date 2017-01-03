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
    double D_e, V, delta, tau_ell, r, a0, a1, p, q, H_bar;
  public:
    MigrationIntegrand(double D_e_, double V_,
                       double delta_, double tau_ell_,
                       double r_,
                       double a0_, double a1_,
                       double p_, double q_, double H_bar_): D_e(D_e_), V(V_),
                                                             delta(delta_),
                                                             tau_ell(tau_ell_),
                                                             r(r_), a0(a0_),
                                                             a1(a1_), p(p_),
                                                             q(q_),
                                                             H_bar(H_bar_) {}

    double operator()(const double& h) const
    {
      double Q = get_Q(V, D_e, delta);
      double upper = (exp(tau_ell * Q - r * h) - 1);
      double log_F_c_given_h = R::pbeta(upper, a0 + a1 * h, 1, 1, 1);
      double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
      return exp(log_F_c_given_h + log_f_h);
    }
};


// [[Rcpp::export]]
double get_Dstar(double D_e, double V, double delta, double tau_ell,
                 double tau_n, double r, double a0, double a1,
                 double p, double q, double H_bar, double omega_n)
{
  // Landholding families
  MigrationIntegrand g(D_e, V, delta, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  const double Q = get_Q(V, D_e, delta);
  const double landed_upper = tau_ell * Q / r;
  const double D_ell = integrate(g, 0.0, landed_upper, err_est, err_code);

  // Landless Families
  const double landless_upper = (exp(tau_n * Q) - 1);
  const double D_n =  R::pbeta(landless_upper, a0, 1, 1, 0);
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
    double D_e, V, delta, tau_ell, r, a0, a1, p, q, H_bar;
  public:
    ExpropriationIntegrand(double D_e_, double V_,
                           double delta_, double tau_ell_,
                           double r_,
                           double a0_, double a1_,
                           double p_, double q_, double H_bar_): D_e(D_e_), V(V_),
                                                                delta(delta_),
                                                                tau_ell(tau_ell_),
                                                                r(r_), a0(a0_),
                                                                a1(a1_), p(p_),
                                                                q(q_),
                                                                H_bar(H_bar_) {}
    double operator()(const double& h) const
    {
      double Q = get_Q(V, D_e, delta);
      double upper = (exp(tau_ell * Q - r * h) - 1);
      double log_F_c_given_h = R::pbeta(upper, a0 + a1 * h, 1, 1, 1);
      double log_f_h = R::dbeta(h / H_bar, p, q, 1) - log(H_bar);
      return exp(log(h) + log_F_c_given_h + log_f_h);
    }
};

// [[Rcpp::export]]
double get_surplus(double V, double delta, double tau_ell, double tau_n,
                   double r, double a0, double a1, double p, double q,
                   double H_bar, double omega_n, double gamma, double beta,
                   double alpha)
{
  double Dstar = get_migration_eq(V, delta, tau_ell, tau_n, r, a0,
                        a1, p, q, H_bar, omega_n);
  ExpropriationIntegrand g(Dstar, V, delta, tau_ell, r, a0, a1, p, q, H_bar);
  double err_est;
  int err_code;
  const double Q = get_Q(V, Dstar, delta);
  const double upper = tau_ell * Q / r;
  const double res = integrate(g, 0.0, upper, err_est, err_code);
  double Xstar = (1 - omega_n) * res;
  return Xstar - gamma * Dstar - beta - alpha * 0.5 * pow(V, 2.0);
}
