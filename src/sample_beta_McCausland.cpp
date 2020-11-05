#include <RcppArmadillo.h>
using namespace Rcpp;
#include "cpp_utilities.h"

void FFBS(arma::mat& beta_nc,
          const arma::vec& y,
          const arma::mat& x,
          const arma::vec& theta_sr,
          const arma::vec& beta_mean,
          const arma::vec& sigma2,
          arma::vec& m_N,
          arma::mat& chol_C_N_inv){

  int d = x.n_cols;
  int N = y.n_elem;

  // Storage for posterior moments
  arma::mat e_store(d, N + 1, arma::fill::zeros);
  arma::cube P_store(d, d, N + 1, arma::fill::zeros);

  // Initial is non-centered so, mean zero & vcov = I
  arma::mat Id = arma::eye(d, d);
  arma::mat Q = Id;
  P_store.slice(0) = Id;

  // Augmented data
  arma::mat x_star = x;
  for (int j = 0; j < d; j++) {
    x_star.col(j) *= theta_sr(j);
  }

  arma::vec y_star = y - x * beta_mean;

  // Forward loop
  for (int t = 1; t < (N + 1); t++) {
    arma::mat P_cond = P_store.slice(t - 1) + Q;
    arma::mat K = P_cond * (x_star.row(t - 1)).t() * arma::as_scalar(1.0/(x_star.row(t - 1) * P_cond * (x_star.row(t - 1)).t() + sigma2(t - 1)));

    P_store.slice(t) = (Id - K * x_star.row(t - 1)) * P_cond;
    e_store.col(t) = e_store.col(t - 1) + K * (y_star(t - 1) - arma::as_scalar(x_star.row(t - 1) * e_store.col(t - 1)));
  }

  // Sample backwards
  arma::vec eps = Rcpp::rnorm(d, 0, 1);
  beta_nc.col(N) = e_store.col(N) + robust_chol(P_store.slice(N)) * eps;


  for (int t = N-1; t >= 0; t--){
    arma::vec mu_post = e_store.col(t) - P_store.slice(t) * arma::solve(P_store.slice(t) + Q, e_store.col(t)) + P_store.slice(t) *
      arma::solve(P_store.slice(t) + Q, beta_nc.col(t + 1));
    arma::mat V_post = P_store.slice(t) - P_store.slice(t) * arma::solve(P_store.slice(t) + Q, P_store.slice(t));

    eps = Rcpp::rnorm(d, 0, 1);
    beta_nc.col(t) = mu_post + robust_chol(V_post) * eps;
  }

  m_N = e_store.col(N);
  arma::mat C_N_inv;
  bool solved = arma::inv_sympd(C_N_inv, P_store.slice(N));

  double jitter = 1e-12 * arma::mean((P_store.slice(N)).diag());
  int max_tries = 1000;
  int num_tries = 1;

  while((solved == false) & (num_tries <= max_tries) & !std::isinf(jitter)) {
    solved = arma::inv_sympd(C_N_inv, P_store.slice(N) + jitter * arma::eye(arma::size(P_store.slice(N))));

    jitter *= 1.1;
    num_tries += 1;
  }
  chol_C_N_inv = robust_chol(C_N_inv);
}

void sample_beta_McCausland(arma::mat& beta_nc_samp,
                            const arma::vec& y,
                            const arma::mat& x,
                            const arma::vec& theta_sr,
                            const arma::vec& sigma2,
                            const arma::vec& beta_mean,
                            arma::vec& m_N,
                            arma::mat& chol_C_N_inv) {

  // This is a pared-down version of the McCausland et al. (2011) algorithm
  // It's fast as hell, but tends to be a bit less stable than classic FFBS
  // Hence it is tried first and supplemented with FFBS in case of failure
  try {
    int d = x.n_cols;
    int N = y.n_elem;

    // helpers
    arma::mat I_d = arma::eye(d, d);

    // Storage objects (for calculating alpha)
    arma::cube L_upper_store(d, d, N+1);
    arma::cube steptwo_store(d, d, N+1);
    arma::cube m_store(d, 1, N+1);

    // Augmented data
    arma::vec y_star = y - x * beta_mean;

    arma::mat x_star = x;
    for (int j = 0; j < d; j++) {
      x_star.col(j) *= theta_sr(j);
    }

    // Omega Mat
    arma::mat Omega = 2*I_d;
    arma::mat Omega_offdiag = -1 * I_d;

    // Vector c
    arma::vec c(d, arma::fill::zeros);

    // Sigma_inverse
    arma::mat Sigma_inv(d, d, arma::fill::zeros);

    // Calculate necessary objects for state 0
    // These are slightly different for the 0th state because some never change for different data sets
    arma::mat L_lower = arma::diagmat(arma::diagvec(Omega)/std::sqrt(2));
    arma::mat L_upper = L_lower.t();

    arma::mat steptwo = -arma::diagmat(arma::diagvec(I_d)/std::sqrt(2));
    arma::mat stepthree = steptwo.t() * steptwo;

    arma::vec a0(d, arma::fill::zeros);
    arma::vec m(d, arma::fill::zeros);

    // Store objects
    L_upper_store.slice(0) = L_upper;
    steptwo_store.slice(0) = steptwo;
    m_store.slice(0) = m;

    // Start forward loop
    for (int t = 1; t < N + 1; t++) {
      Omega = (x_star.row(t-1)).t() * x_star.row(t-1) * 1.0/sigma2(t-1) + (1 + (t != N))*I_d;
      c = (x_star.row(t-1)).t()/sigma2(t-1) * y_star(t-1);
      Sigma_inv = Omega - stepthree;

      L_lower = robust_chol(Sigma_inv);
      L_upper = L_lower.t();

      steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);

      stepthree = steptwo.t() * steptwo;

      m = arma::solve(arma::trimatu(L_upper), arma::solve(arma::trimatl(L_lower), (c - Omega_offdiag * m)));
      L_upper_store.slice(t) = L_upper;
      steptwo_store.slice(t) = steptwo;
      m_store.slice(t) = m;
    }


    arma::vec eps = Rcpp::rnorm(d, 0, 1);
    arma::mat l = arma::solve(arma::trimatu(L_upper_store.slice(N)), eps);
    beta_nc_samp.col(N) = l + m_store.slice(N);

    for (int t = N-1; t >= 0; t--){
      eps = Rcpp::rnorm(d, 0, 1);
      arma::vec q = eps - steptwo_store.slice(t) * beta_nc_samp.col(t+1);
      l = arma::solve(arma::trimatu(L_upper_store.slice(t)), q);
      beta_nc_samp.col(t) = l + m_store.slice(t);
    }

    std::for_each(beta_nc_samp.begin(), beta_nc_samp.end(), res_protector);


    // Return objects for LPDS calculation
    m_N = m_store.slice(N);
    chol_C_N_inv = (L_upper_store.slice(N)).t();
  } catch (...) {
    // Fall back on FFBS if McCausland fails
    FFBS(beta_nc_samp,
         y,
         x,
         theta_sr,
         beta_mean,
         sigma2,
         m_N,
         chol_C_N_inv);
  }

}
