// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppNumerical)]]
#include <RcppNumerical.h>
using namespace Numer;

class MigrationIntegrand: public Func
{
  private:
    double D_e;
    double r;
    double V;
    double P;
    double s_h;
    double mu_h;
    double tau_ell;
    double s_c;
    double mu_c;
  public:
    MigrationIntegrand(double D_e_, double r_, double V_, double P_,
                       double s_h_, double mu_h_, double tau_ell_,
                       double s_c_, double mu_c_) : D_e(D_e_),  r(r_), V(V_),
                                                    P(P_),  s_h(s_h_),
                                                    mu_h(mu_h_),
                                                    tau_ell(tau_ell_),
                                                    s_c(s_c_), mu_c(mu_c_) {}
    double operator()(const double& h) const
    {
      const double Q = V / (P * (1 - D_e));
      double F_c = R::plnorm(tau_ell * Q - r * h, mu_c, s_c, 1, 0);
      double f_c = R::dlnorm(h, mu_h, s_h, 0);
      return F_c * f_c;
    }

};

// [[Rcpp::export]]
double get_Dstar(double D_e, double r, double V, double P, double s_h, double mu_h,
             double tau_ell, double s_c, double mu_c, double tau_n,
             double delta, double frac_n)
{
  // Landholding families
  MigrationIntegrand g(D_e, r, V, P, s_h, mu_h, tau_ell, s_c, mu_c);
  double err_est;
  int err_code;
  double upper = tau_ell * (V / (P * (1 - D_e))) / r;
  const double res = integrate(g, 0.0, upper, err_est, err_code);
  const double D_ell = (1 - delta) * res;

  // Landless Families
  const double Q = V / (P * (1 - D_e));
  const double D_n = (1 - delta) * R::plnorm(tau_n * Q, mu_c, s_c, 1, 0);

  return frac_n * D_n + (1 - frac_n) * D_ell;

}

// [[Rcpp::export]]
double get_migration_eq(double r, double V, double P, double s_h, double mu_h,
                  double tau_ell, double s_c, double mu_c, double tau_n,
                  double delta, double frac_n){
  const int max_iter = 500;
  const double tol = 0.0001;
  int i = 0;
  double x = 0.0;
  double f = get_Dstar(x, r, V, P, s_h, mu_h, tau_ell, s_c, mu_c, tau_n, delta, frac_n);
  while((i < max_iter) & (std::abs(f - x) > tol)){
    x = f;
    f = get_Dstar(x, r, V, P, s_h, mu_h, tau_ell, s_c, mu_c, tau_n, delta, frac_n);
    i++;
  }
  return(x);
}

