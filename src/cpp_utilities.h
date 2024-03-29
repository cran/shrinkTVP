#ifndef CPP_UTILITIES_H
#define CPP_UTILITIES_H

void res_protector(double& x);

void calc_xi2_tau2(arma::vec& param,
                   const arma::vec& param_til,
                   const arma::vec& loc_shrink_til,
                   double glob_shrink,
                   double c,
                   double a);

void to_CP(arma::mat& beta,
           const arma::mat& beta_nc,
           const arma::vec& theta_sr,
           const arma::vec& beta_mean);

void to_NCP(arma::mat& beta_nc,
            const arma::mat& beta,
            const arma::vec& theta_sr,
            const arma::vec& beta_mean);

arma::mat robust_chol (const arma::mat& V);

arma::mat robust_chol_nontri (const arma::mat& V);

arma::mat robust_solve(arma::mat A,
                       arma::mat B);

double unur_bessel_k_nuasympt(double x,
                              double nu,
                              bool islog,
                              bool expon_scaled);

void sample_lin_reg_stab(arma::vec& param_vec,
                         const arma::vec& y,
                         const arma::mat& X,
                         const arma::vec& sigma2,
                         const arma::vec& prior_var);

void sample_lin_reg_rue(arma::vec& param_vec,
                        const arma::vec& y,
                        const arma::mat& X,
                        const arma::vec& sigma2,
                        const arma::vec& prior_var);

void sample_lin_reg_rue_homosc(arma::vec& param_vec,
                               const arma::vec& Xty,
                               const arma::mat& XtX,
                               double sigma2,
                               const arma::vec& prior_var);

void sample_lin_reg_bhat(arma::vec& param_vec,
                         const arma::vec& y,
                         const arma::mat& x,
                         double sigma2,
                         const arma::vec& prior_var);

double samp_disc_given(arma::vec to_sample,
                       arma::vec probs);

double samp_disc(const arma::rowvec& lambda_p,
                 double a,
                 double c,
                 double rho,
                 double tol,
                 int max_size);

double samp_disc2(double kappa_pre,
                  double kappa_post,
                  double psi_curr,
                  double psi_pre,
                  double a,
                  double c,
                  double rho,
                  double tol,
                  int max_size);

#endif
