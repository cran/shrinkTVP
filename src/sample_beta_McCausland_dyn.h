#ifndef SAMPLE_BETA_MCCAUSLAND_DYN
#define SAMPLE_BETA_MCCAUSLAND_DYN

void sample_beta_McCausland_dyn(arma::mat& beta_nc_samp,
                                const arma::vec& y,
                                const arma::mat& x,
                                const arma::vec& theta_sr,
                                const arma::vec& sigma2,
                                const arma::vec& beta_mean,
                                const arma::mat& psi,
                                arma::vec& m_N,
                                arma::mat& chol_C_N_inv);

#endif
