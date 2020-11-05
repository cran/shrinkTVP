#ifndef SAMPLE_BETA_MCCAUSLAND_H
#define SAMPLE_BETA_MCCAUSLAND_H

void FFBS(arma::mat& beta_nc,
          const arma::vec& y,
          const arma::mat& x,
          const arma::vec& theta_sr,
          const arma::vec& beta_mean,
          const arma::vec& sigma2,
          arma::vec& m_N,
          arma::mat& chol_C_N_inv);

void sample_beta_McCausland(arma::mat& beta_nc_samp,
                                const arma::vec& y,
                                const arma::mat& x,
                                const arma::vec& theta_sr,
                                const arma::vec& sigma2,
                                const arma::vec& beta_mean,
                                arma::vec& m_N,
                                arma::mat& chol_C_N_inv);

#endif
