#ifndef COMMON_SAMPLING_FUNCTIONS_H
#define COMMON_SAMPLING_FUNCTIONS_H

void sample_alpha(arma::vec& beta_mean,
                  arma::vec& theta_sr,
                  const arma::vec& y,
                  const arma::mat& x,
                  const arma::mat& beta_nc,
                  const arma::vec& sigma2,
                  const arma::vec& tau2,
                  const arma::vec& xi2);

void resample_alpha(arma::vec& beta_mean,
                    arma::vec& theta_sr,
                    const arma::mat& beta,
                    const arma::mat& beta_nc,
                    const arma::vec& xi2,
                    const arma::vec& tau2);

void sample_sigma2(arma::vec& sigma2,
                   const arma::vec& y,
                   const arma::mat& x,
                   const arma::mat beta_nc,
                   const arma::vec beta_mean,
                   const arma::vec& theta_sr,
                   double c0,
                   double C0);

double sample_C0(const arma::vec& sigma2,
                 double g0,
                 double c0,
                 double G0);

#endif
