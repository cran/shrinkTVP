// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "do_rgig1.h"
#include <math.h>
using namespace Rcpp;


void res_protector(double& x){
  if (std::abs(x) < DBL_MIN * std::pow(10, 10)){
    double sign = std::copysign(1, x);
    x = DBL_MIN * std::pow(10, 10) * sign;
  }
}


void sample_alpha(arma::vec& alpha_samp, arma::vec& y, arma::mat& x, arma::mat& W, arma::colvec& tau2, arma::colvec& xi2, arma::vec& sig2, arma::vec& a0, int d, Function Rchol) {
  arma::mat W_til = W.t() * arma::diagmat(1/sig2);
  arma::mat A0_sr = arma::diagmat(arma::sqrt(arma::join_cols(tau2, xi2)));
  arma::vec prior_a = a0/arma::join_cols(tau2, xi2);
  arma::mat a = W_til * y + prior_a;
  arma::mat Omega_star = A0_sr * W_til * W * A0_sr + arma::eye(2*d, 2*d);

  arma::mat A_t;
  arma::mat A_t_til;
  bool solved = arma::solve(A_t_til, Omega_star, A0_sr, arma::solve_opts::no_approx);
  if (solved == 1){
    A_t = A0_sr * A_t_til;
  } else {
    arma::mat A0 = arma::diagmat(arma::join_cols(tau2, xi2));
    A_t = arma::inv(W_til * W + arma::inv(diagmat(A0)));
  }

  arma::vec v = rnorm(2*d);

  arma::mat cholA;
  bool chol_success = chol(cholA, A_t);

  // Fall back on Rs chol if armadillo fails (it suppports pivoting)
  if (chol_success == false){
    Rcpp::NumericMatrix tmp = Rchol(A_t, true, false, -1);
    arma::mat cholA_tmp = arma::mat(tmp.begin(), 2*d, 2*d, false);
    arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
    cholA = cholA_tmp.cols(piv);
  }

  alpha_samp = A_t * a + cholA.t() * v;
  std::for_each(alpha_samp.begin(), alpha_samp.end(), res_protector);
}


void resample_alpha_diff(arma::mat& alpha_samp, arma::mat& betaenter, arma::vec& theta_sr, arma::vec& beta_mean, arma::mat& beta_diff,  arma::vec& xi2, arma::vec& tau2, int d, int N){
  arma::vec sign_sqrt = arma::sign(theta_sr);
  arma::colvec theta_sr_new(d, arma::fill::none);
  int p1_theta = -N/2;

  arma::vec theta(d, arma::fill::none);


  for (int j = 0; j < d; j++){
    double p2_theta = 1/xi2(j);
    double p3_theta = arma::as_scalar(arma::accu(arma::pow(beta_diff.row(j), 2))) +
      std::pow((betaenter(j, 0) - beta_mean(j)), 2);

    double res = do_rgig1(p1_theta, p3_theta, p2_theta);
    theta(j) = res;
    theta_sr_new(j) = std::sqrt(arma::as_scalar(theta(j))) * sign_sqrt(j);
  }

  arma::colvec beta_mean_new(d, arma::fill::none);

  for (int j = 0; j < d; j++){
    double sigma2_beta_mean = 1/(1/tau2(j) + 1/(theta(j)));
    double mu_beta_mean = betaenter(j, 0) * tau2(j)/(tau2(j) + theta(j));
    beta_mean_new(j) = R::rnorm(mu_beta_mean, std::sqrt(sigma2_beta_mean));
  }

  alpha_samp = arma::join_cols(beta_mean_new, theta_sr_new);
  std::for_each(alpha_samp.begin(), alpha_samp.end(), res_protector);
}




void sample_xi2(arma::vec& xi2_samp, arma::vec& theta_sr, double kappa2, double a_xi, int d){
  arma::vec theta = arma::pow(theta_sr, 2);
  arma::vec xi2(d, arma::fill::none);

  for (int j = 0; j < d; j++){

    double p1_xi = a_xi - 0.5;
    double p2_xi = a_xi * kappa2;
    double p3_xi = theta(j);

    double res = do_rgig1(p1_xi, p3_xi, p2_xi);

    res_protector(res);

    xi2_samp(j) = res;
  }

}


void sample_tau2(arma::vec& tau2_samp, arma::vec& beta_mean, double lambda2, double a_tau, int d){
  arma::vec tau2(d, arma::fill::none);

  for (int j = 0; j < d; j++){
    double p1_tau = a_tau - 0.5;
    double p2_tau = a_tau * lambda2;
    double p3_tau = std::pow(beta_mean(j), 2);

    double res = do_rgig1(p1_tau, p3_tau, p2_tau);

    res_protector(res);

    tau2_samp(j) = res;

  }
}


double sample_kappa2(arma::vec& xi2, double a_xi, double d1, double d2, int d){
  double d1_full = d1 + a_xi * d;
  double d2_full = d2 + arma::as_scalar(arma::mean(xi2)) * a_xi * d * 0.5;
  double kappa2 = R::rgamma(d1_full, 1/d2_full);
  return(kappa2);
}


double sample_lambda2(arma::vec& tau2, double a_tau, double e1, double e2, int d){
  double e1_full = e1 + a_tau * d;
  double e2_full = e2 + arma::as_scalar(arma::mean(tau2)) * a_tau * d * 0.5;
  double lambda2 = R::rgamma(e1_full, 1/e2_full);
  return(lambda2);
}


void sample_sigma2(arma::vec& sig2_samp, arma::vec& y, arma::mat& W, arma::vec& alpha, double c0, double C0, int N){
  double a_full = c0 + N/2;
  double b_full = C0 + 0.5 * arma::as_scalar(arma::sum(arma::pow((y - W*alpha), 2)));
  double sig2 = 1/R::rgamma(a_full, 1/b_full);
  sig2_samp.fill(sig2);
}


double sample_C0(arma::vec& sig2, double g0, double c0, double G0){
  double a_full = g0 + c0;
  double b_full = G0 + 1/sig2(0);
  double C0 = R::rgamma(a_full, 1/b_full);
  return(C0);
}



