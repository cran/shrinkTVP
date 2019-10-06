// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;


void sample_beta_McCausland(arma::mat& beta_nc_samp, arma::vec& y, arma::mat& x,
                            arma::colvec& theta_sr, arma::vec& sig2, arma::colvec& beta_mean,
                            arma::vec& m_N, arma::mat& chol_C_N_inv, bool LPDS, int N, int d,
                            Function Rchol) {

  arma::colvec theta = arma::pow(theta_sr, 2);

  arma::cube Omega_diag(d, d, N+1);
  arma::mat I_d = arma::eye(d, d);
  arma::mat theta_sr_diag = arma::diagmat(theta_sr);

  Omega_diag.slice(0) = 2 * I_d;
  for (int t = 1; t < N; t++){
    Omega_diag.slice(t) = (x.row(t-1) * (theta_sr_diag)).t() * (x.row(t-1)*(theta_sr_diag % I_d))/sig2(t-1) + 2*I_d;
  }
  Omega_diag.slice(N) = (x.row(N-1) * (theta_sr_diag)).t() * (x.row(N-1)*(theta_sr_diag))/sig2(N-1) + I_d;

  arma::mat Omega_offdiag = -1 * I_d;
  arma::vec y_star = y - x * beta_mean;
  arma::cube cvett(d, 1, N+1);
  cvett.slice(0) = arma::vec(d, arma::fill::zeros);

  for (int t = 1; t < N+1; t++){
    cvett.slice(t) = (x.row(t-1) * (theta_sr_diag * I_d)).t() / sig2(t-1) * y_star(t-1);
  }

  arma::mat Sigma_inv = Omega_diag.slice(0);

  arma::mat L_upper;
  bool chol_success = chol(L_upper, Sigma_inv);

  // Fall back on Rs chol if armadillo fails (it suppports pivoting)
  if (chol_success == false){
    Rcpp::NumericMatrix tmp = Rchol(Sigma_inv, true, false, -1);
    arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
    arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
    L_upper = L_upper_tmp.cols(piv);
  }

  arma::mat L_lower = L_upper.t();

  arma::cube L_upper_list(d, d, N+1);
  arma::cube L_lower_list(d, d, N+1);

  L_upper_list.slice(0) = L_upper;
  L_lower_list.slice(0) = L_lower;

  arma::mat steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);
  arma::cube steptwo_list(d, d, N+1);
  steptwo_list.slice(0) = steptwo;

  arma::mat stepthree = steptwo.t() * steptwo;

  arma::mat a0 = solve(arma::trimatl(L_lower), cvett.slice(0));
  arma::mat m = solve(arma::trimatu(L_upper), a0);
  arma::cube m_list(d, 1, N+1);
  m_list.slice(0) = m;

  for (int t = 1; t < N+1; t++){
    Sigma_inv = Omega_diag.slice(t) - stepthree;
    chol_success = chol(L_upper, Sigma_inv);

    // Fall back on Rs chol if armadillo fails (it suppports pivoting)
    if (chol_success == false){
      Rcpp::NumericMatrix tmp = Rchol(Sigma_inv, true, false, -1);
      arma::uvec piv = arma::sort_index(as<arma::vec>(tmp.attr("pivot")));
      arma::mat L_upper_tmp = arma::mat(tmp.begin(), d, d, false);
      L_upper = L_upper_tmp.cols(piv);
    }
    L_lower = L_upper.t();

    L_upper_list.slice(t) = L_upper;
    L_lower_list.slice(t) = L_lower;

    steptwo = solve(arma::trimatl(L_lower), Omega_offdiag);
    steptwo_list.slice(t) = steptwo;

    stepthree = steptwo.t() * steptwo;

    arma::mat l = cvett.slice(t) - Omega_offdiag.t() * m_list.slice(t-1);
    arma::mat at = arma::solve(arma::trimatl(L_lower), l);
    m = arma::solve(arma::trimatu(L_upper), at);
    m_list.slice(t) = m;
  }

  arma::vec eps = Rcpp::rnorm(d, 0, 1);
  arma::mat l = arma::solve(arma::trimatu(L_upper_list.slice(N)), eps);
  beta_nc_samp.col(N) = l + m_list.slice(N);

  for (int t = N-1; t >= 0; t--){
    eps = Rcpp::rnorm(d, 0, 1);
    arma::vec q = eps - steptwo_list.slice(t) * beta_nc_samp.col(t+1);
    l = arma::solve(arma::trimatu(L_upper_list.slice(t)), q);
    beta_nc_samp.col(t) = l + m_list.slice(t);
  }

  if (LPDS == true){
    m_N = m_list.slice(N);
    chol_C_N_inv = L_lower_list.slice(N);
  }

}
