#include <RcppArmadillo.h>
#include "do_rgig1.h"
#include "cpp_utilities.h"
#include <math.h>
using namespace Rcpp;



void sample_alpha(arma::vec& beta_mean,
                  arma::vec& theta_sr,
                  const arma::vec& y,
                  const arma::mat& x,
                  const arma::mat& beta_nc,
                  const arma::vec& sigma2,
                  const arma::vec& tau2,
                  const arma::vec& xi2) {

  int d = x.n_cols;
  int N = y.n_elem;
  arma::mat x_tilde = x % (beta_nc.cols(1,N)).t();
  arma::mat W = arma::join_rows(x, x_tilde);

  arma::vec alpha_samp(2*d);

  sample_lin_reg_stab(alpha_samp,
                      y,
                      W,
                      sigma2,
                      arma::join_cols(tau2, xi2));


  std::for_each(alpha_samp.begin(), alpha_samp.end(), res_protector);

  beta_mean = alpha_samp.rows(0, d-1);
  theta_sr = alpha_samp.rows(d, 2*d-1);
}


void resample_alpha(arma::vec& beta_mean,
                    arma::vec& theta_sr,
                    const arma::mat& beta,
                    const arma::mat& beta_nc,
                    const arma::vec& xi2,
                    const arma::vec& tau2){

  // Difference beta beforehand (for numerical stability)
  arma::mat beta_diff_pre = arma::diff(beta_nc, 1, 1);
  arma::mat beta_diff =  beta_diff_pre.each_col() % theta_sr;

  // Get some necessary dimensions and set up storage items
  int d = theta_sr.n_elem;
  int N = beta_nc.n_cols - 1;
  arma::vec sign_sqrt = arma::sign(theta_sr);
  arma::vec theta(d, arma::fill::none);
  arma::colvec theta_sr_new(d, arma::fill::none);
  arma::colvec beta_mean_new(d, arma::fill::none);

  // The first parameter is the same across all thetas
  int p1_theta = -N * 0.5;

  // Sample theta in centered parameterization
  for (int j = 0; j < d; j++){
    double p2_theta = 1.0/xi2(j);
    double p3_theta = arma::as_scalar(arma::accu(arma::pow(beta_diff.row(j), 2))) +
      std::pow((beta(j, 0) - beta_mean(j)), 2);

    theta(j) = do_rgig1(p1_theta, p3_theta, p2_theta);
    theta_sr_new(j) = std::sqrt(theta(j)) * sign_sqrt(j);
  }

  // Sample beta in centered parameterization
  for (int j = 0; j < d; j++){
    double sigma2_beta_mean = 1.0/(1.0/tau2(j) + 1.0/(theta(j)));
    double mu_beta_mean = beta(j, 0) * tau2(j)/(tau2(j) + theta(j));
    beta_mean_new(j) = R::rnorm(mu_beta_mean, std::sqrt(sigma2_beta_mean));
  }

  std::for_each(theta_sr_new.begin(), theta_sr_new.end(), res_protector);
  std::for_each(beta_mean_new.begin(), beta_mean_new.end(), res_protector);

  beta_mean = beta_mean_new;
  theta_sr = theta_sr_new;
}


void sample_sigma2(arma::vec& sigma2,
                   const arma::vec& y,
                   const arma::mat& x,
                   const arma::mat beta_nc,
                   const arma::vec beta_mean,
                   const arma::vec& theta_sr,
                   double c0,
                   double C0){

  int N = y.n_elem;
  arma::vec alpha = arma::join_cols(beta_mean, theta_sr);
  arma::mat x_tilde = x % (beta_nc.cols(1,N)).t();
  arma::mat W = arma::join_rows(x, x_tilde);

  double a_full = c0 + N * 0.5;
  double b_full = C0 + 0.5 * arma::as_scalar(arma::sum(arma::pow((y - W*alpha), 2)));
  double sigma2_double = 1.0/R::rgamma(a_full, 1.0/b_full);

  res_protector(sigma2_double);
  sigma2.fill(sigma2_double);
}


double sample_C0(const arma::vec& sigma2,
                 double g0,
                 double c0,
                 double G0){

  double a_full = g0 + c0;
  double b_full = G0 + 1.0/sigma2(0);
  double C0 = R::rgamma(a_full, 1.0/b_full);
  return(C0);
}
