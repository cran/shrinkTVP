#ifndef DYN_SAMPLING_STEPS_H
#define DYN_SAMPLING_STEPS_H

double sample_lambda_iid(double a,
                         double c,
                         double psi);

double sample_psi_iid(double c,
                      double lambda,
                      double w);

double sample_lambda_0(double kappa_1,
                       double a,
                       double c,
                       double rho);

arma::rowvec sample_lambda(const arma::vec & kappa,
                           const arma::vec & psi,
                           double a,
                           double c,
                           double rho);

arma::rowvec sample_psi(const arma::vec& lambda,
                        const arma::vec& beta_nc,
                        double c);

void sample_kappa_fast_marg_alternating(arma::vec& kappa,
                                        const arma::vec& psi,
                                        double a,
                                        double c,
                                        double rho);
#endif
