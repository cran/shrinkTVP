#include <RcppArmadillo.h>
using namespace Rcpp;
#include "cpp_utilities.h"

void FFBS(arma::mat& beta_nc,
          const arma::vec& y,
          const arma::mat& x,
          const arma::vec& theta_sr,
          const arma::vec& beta_mean,
          const arma::vec& sigma2,
          const arma::mat& psi,
          arma::vec& m_N,
          arma::mat& chol_C_N_inv){

  int d = x.n_cols;
  int N = y.n_elem;

  // Storage for posterior moments
  arma::mat e_store(d, N, arma::fill::zeros);
  arma::cube P_store(d, d, N, arma::fill::zeros);

  // Augmented data
  arma::mat x_star = x;
  for (int j = 0; j < d; j++) {
    x_star.col(j) *= theta_sr(j);
  }

  arma::vec y_star = y - x * beta_mean;


  // Initial is non-centered so, mean zero & vcov = psi(0)
  arma::mat Id = arma::eye(d, d);
  arma::mat Q = Id;

  Q = arma::diagmat(psi.col(0));
  arma::mat P_cond = Q;
  arma::mat K = P_cond * (x_star.row(0)).t() * arma::as_scalar(1.0/(x_star.row(0) * P_cond * (x_star.row(0)).t() + sigma2(0)));

  P_store.slice(0) = (Id - K * x_star.row(0)) * P_cond;
  e_store.col(0) = K * y_star(0);


  // Forward loop
  for (int t = 1; t < N; t++) {
    Q = arma::diagmat(psi.col(t));
    arma::mat P_cond = P_store.slice(t - 1) + Q;
    arma::mat K = P_cond * (x_star.row(t)).t() * arma::as_scalar(1.0/(x_star.row(t) * P_cond * (x_star.row(t)).t() + sigma2(t)));

    P_store.slice(t) = (Id - K * x_star.row(t)) * P_cond;
    e_store.col(t) = e_store.col(t - 1) + K * (y_star(t) - arma::as_scalar(x_star.row(t) * e_store.col(t - 1)));
  }

  // Sample backwards
  arma::vec eps = Rcpp::rnorm(d, 0, 1);
  beta_nc.col(N-1) = e_store.col(N-1) + robust_chol_nontri(P_store.slice(N-1)) * eps;


  for (int t = N-2; t >= 0; t--){
    Q = arma::diagmat(psi.col(t+1));
    arma::vec mu_post = e_store.col(t) - P_store.slice(t) * arma::solve(P_store.slice(t) + Q, e_store.col(t)) +
      P_store.slice(t) * arma::solve(P_store.slice(t) + Q, beta_nc.col(t + 1));
    arma::mat V_post = P_store.slice(t) - P_store.slice(t) * arma::solve(P_store.slice(t) + Q, P_store.slice(t));

    eps = Rcpp::rnorm(d, 0, 1);
    beta_nc.col(t) = mu_post + robust_chol_nontri(V_post) * eps;
  }

  m_N = e_store.col(N-1);
  arma::mat C_N_inv;
  bool solved = arma::inv_sympd(C_N_inv, P_store.slice(N-1));

  double jitter = 1e-12 * arma::mean((P_store.slice(N-1)).diag());
  int max_tries = 1000;
  int num_tries = 1;

  while((solved == false) & (num_tries <= max_tries) & !std::isinf(jitter)) {
    solved = arma::inv_sympd(C_N_inv, P_store.slice(N-1) + jitter * arma::eye(arma::size(P_store.slice(N-1))));

    jitter *= 1.1;
    num_tries += 1;
  }
  chol_C_N_inv = robust_chol(C_N_inv);
}

void sample_beta_McCausland_dyn(arma::mat& beta_nc_samp,
                                const arma::vec& y,
                                const arma::mat& x,
                                const arma::vec& theta_sr,
                                const arma::vec& sigma2,
                                const arma::vec& beta_mean,
                                const arma::mat& psi,
                                arma::vec& m_N,
                                arma::mat& chol_C_N_inv) {

  try {
    // This is a pared-down version of the McCausland et al. (2011) algorithm
    int d = x.n_cols;
    int N = y.n_elem;

    // helpers
    arma::mat I_d = arma::eye(d, d);
    arma::mat theta_sr_diag = arma::diagmat(theta_sr);

    // Storage objects (for calculating alpha)
    arma::cube L_upper_store(d, d, N);
    arma::cube steptwo_store(d, d, N);
    arma::cube m_store(d, 1, N);

    // Augmented data
    arma::vec y_star = y - x * beta_mean;
    arma::mat x_star = x * theta_sr_diag;

    // Omega Mats
    arma::mat psi_inv = 1/psi;
    arma::mat Omega = (x_star.row(0)).t() * x_star.row(0) * 1/sigma2(0) +
      arma::diagmat(psi_inv.col(0)) + arma::diagmat(psi_inv.col(1));
    arma::mat Omega_pre = -I_d;
    arma::mat Omega_post = -arma::diagmat(psi_inv.col(1));
    arma::mat Omega_final;

    // Vector c
    arma::vec c = (x_star.row(0)).t()/sigma2(0) * y_star(0);

    // Sigma_inverse
    arma::mat Sigma_inv(d, d, arma::fill::zeros);

    // Calculate necessary objects for state 0
    arma::mat L_lower = arma::chol(Omega, "lower");
    arma::mat L_upper = L_lower.t();

    arma::mat steptwo = solve(arma::trimatl(L_lower), Omega_post);
    arma::mat stepthree = steptwo.t() * steptwo;

    arma::mat a0 = solve(arma::trimatl(L_lower), c);
    arma::mat m = solve(arma::trimatu(L_upper), a0);

    // Store objects
    L_upper_store.slice(0) = L_upper;
    steptwo_store.slice(0) = steptwo;
    m_store.slice(0) = m;


    // Start forward loop
    for (int t = 1; t < N; t++) {
      Omega_pre = -arma::diagmat(psi_inv.col(t));

      if (t != (N-1)) {
        Omega_final = arma::diagmat(psi_inv.col(t+1));
        Omega_post = -arma::diagmat(psi_inv.col(t+1));
      } else {
        Omega_final.fill(0);
      }

      Omega = (x_star.row(t)).t() * x_star.row(t) * 1/sigma2(t) + arma::diagmat(psi_inv.col(t)) +  Omega_final;
      c = (x_star.row(t)).t()/sigma2(t) * y_star(t);
      Sigma_inv = Omega - stepthree;

      L_lower = arma::chol(Sigma_inv, "lower");
      L_upper = L_lower.t();

      if (t != (N-1)) {
        steptwo = solve(arma::trimatl(L_lower), Omega_post);
        stepthree = steptwo.t() * steptwo;
      }

      m = arma::solve(arma::trimatu(L_upper), arma::solve(arma::trimatl(L_lower), (c - Omega_pre * m)));
      L_upper_store.slice(t) = L_upper;
      steptwo_store.slice(t) = steptwo;
      m_store.slice(t) = m;

    }


    arma::vec eps = Rcpp::rnorm(d, 0, 1);
    arma::mat l = arma::solve(arma::trimatu(L_upper_store.slice(N-1)), eps);
    beta_nc_samp.col(N-1) = l + m_store.slice(N-1);

    for (int t = N-2; t >= 0; t--){
      eps = Rcpp::rnorm(d, 0, 1);
      arma::vec q = eps - steptwo_store.slice(t) * beta_nc_samp.col(t+1);
      l = arma::solve(arma::trimatu(L_upper_store.slice(t)), q);
      beta_nc_samp.col(t) = l + m_store.slice(t);
    }



    // Return objects for LPDS calculation
    m_N = m_store.slice(N-1);
    chol_C_N_inv = (L_upper_store.slice(N-1)).t();
  } catch (...) {
    FFBS(beta_nc_samp,
         y,
         x,
         theta_sr,
         beta_mean,
         sigma2,
         psi,
         m_N,
         chol_C_N_inv);
  }

  std::for_each(beta_nc_samp.begin(), beta_nc_samp.end(), res_protector);
}
