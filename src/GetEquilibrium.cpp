// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]

#include <RcppNumerical.h>
using namespace Numer;

double get_Q(double V, double Gamma){
  return((Gamma / (1 - Gamma)) * (1 / (Gamma + (1 - Gamma) * exp(-V)) - 1));
}

class MigrationIntegrand: public Func
{
  private:
    double D_e, V, m0, m1, tau_ell, r, u_bar, a0, a1, p, q, H;
  public:
    MigrationIntegrand(double D_e_, double V_,
                       double m0_, double m1_,
                       double tau_ell_, double r_,
                       double u_bar_, double a0_,
                       double a1_, double p_,
                       double q_, double H_): D_e(D_e_), V(V_),
                                              m0(m0_), m1(m1_),
                                              tau_ell(tau_ell_), r(r_),
                                              u_bar(u_bar_), a0(a0_),
                                              a1(a1_), p(p_),
                                              q(q_), H(H_) {}

    double operator()(const double& h) const
    {
      double Q = get_Q(V, m0 + m1 * D_e);
      double upper = (exp(tau_ell * Q - r * h) - 1) * exp(-u_bar);
      double F_c_given_h = R::pbeta(upper, a0 + a1 * h, 1, 1, 0);
      double f_h = pow(h, p - 1) * pow(H - h, q - 1) / (R::beta(p, q) * pow(H, p + q));
      return F_c_given_h * f_h;
    }
};


// [[Rcpp::export]]
double get_Dstar(double D_e, double V, double m0, double m1, double tau_ell,
                 double tau_n, double r, double u_bar, double a0, double a1,
                 double p, double q, double H, double delta, double omega_n)
{
  // Landholding families
  MigrationIntegrand g(D_e, V, m0, m1, tau_ell, r, u_bar, a0, a1, p, q, H);
  double err_est;
  int err_code;
  const double Q = get_Q(V, m0 + m1 * D_e);
  const double landed_upper = tau_ell * Q / r;
  const double res = integrate(g, 0.0, landed_upper, err_est, err_code);
  const double D_ell = (1 - delta) * res;

  // Landless Families
  const double landless_upper = (exp(tau_n * Q) - 1) * exp(-u_bar);
  const double D_n = (1 - delta) * R::pbeta(landless_upper, a0, 1, 1, 0);
  return omega_n * D_n + (1 - omega_n) * D_ell;
}


// [[Rcpp::export]]
double get_migration_eq(double V, double m0, double m1, double tau_ell,
                        double tau_n, double r, double u_bar, double a0,
                        double a1, double p, double q, double H, double delta,
                        double omega_n)
{
  const int max_iter = 500;
  const double tol = 0.0001;
  int i = 0;
  double x = 0.0;
  double f = get_Dstar(x, V, m0, m1, tau_ell, tau_n, r, u_bar, a0, a1, p, q, H,
                  delta, omega_n);
  while((i < max_iter) & (std::abs(f - x) > tol)){
    x = f;
    f = get_Dstar(x, V, m0, m1, tau_ell, tau_n, r, u_bar, a0, a1, p, q, H,
                  delta, omega_n);
    i++;
  }
  return(x);
}

class ExpropriationIntegrand: public Func
{
  private:
    double D_e, V, m0, m1, tau_ell, r, u_bar, a0, a1, p, q, H;
  public:
    ExpropriationIntegrand(double D_e_, double V_,
                           double m0_, double m1_,
                           double tau_ell_, double r_,
                           double u_bar_, double a0_,
                           double a1_, double p_,
                           double q_, double H_): D_e(D_e_), V(V_),
                                                  m0(m0_), m1(m1_),
                                                  tau_ell(tau_ell_), r(r_),
                                                  u_bar(u_bar_), a0(a0_),
                                                  a1(a1_), p(p_),
                                                  q(q_), H(H_) {}

    double operator()(const double& h) const
    {
      double Q = get_Q(V, m0 + m1 * D_e);
      double upper = (exp(tau_ell * Q - r * h) - 1) * exp(-u_bar);
      double F_c_given_h = R::pbeta(upper, a0 + a1 * h, 1, 1, 0);
      double f_h = pow(h, p - 1) * pow(H - h, q - 1) / (R::beta(p, q) * pow(H, p + q));
      return h * F_c_given_h * f_h;
    }
};

// [[Rcpp::export]]
double get_surplus(double V, double m0, double m1, double tau_ell, double tau_n,
                   double r, double u_bar, double a0, double a1, double p,
                   double q, double H, double delta, double omega_n,
                   double gamma, double beta)
{
  double Dstar = get_migration_eq(V, m0, m1, tau_ell, tau_n, r, u_bar, a0,
                        a1, p, q, H, delta, omega_n);
  ExpropriationIntegrand g(Dstar, V, m0, m1, tau_ell, r, u_bar, a0, a1, p, q, H);
  double err_est;
  int err_code;
  const double Q = get_Q(V, m0 + m1 * Dstar);
  const double upper = tau_ell * Q / r;
  const double res = integrate(g, 0.0, upper, err_est, err_code);
  double Xstar = (1 - delta) * (1 - omega_n) * res;
  return Xstar - gamma * Dstar - beta;
}
