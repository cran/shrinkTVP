#ifndef SAMPLE_BETA_MCCAUSLAND_H
#define SAMPLE_BETA_MCCAUSLAND_H

#include <RcppArmadillo.h>

using namespace Rcpp;

void sample_beta_McCausland(arma::mat& beta_nc_samp, arma::vec& y, arma::mat& x,
                            arma::colvec& theta_sr, arma::vec& sig2, arma::colvec& beta_mean,
                            arma::vec& m_N, arma::mat& C_N_inv, bool LPDS, int N, int d,
                            Function Rchol);
#endif
