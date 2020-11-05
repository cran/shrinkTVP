#ifndef TG_SAMPLING_FUNCTIONS_H
#define TG_SAMPLING_FUNCTIONS_H


void TG_sample_prior_var_til(arma::vec& prior_var_til,
                             const arma::vec& param_vec,
                             const arma::vec& local_shrink,
                             double global_shrink,
                             double a,
                             double c);

double TG_sample_global_shrink(const arma::vec& prior_var_til,
                               const arma::vec& local_shrink,
                               const arma::vec& param_vec,
                               double a,
                               double c,
                               double hyper2,
                               bool common_shrink_par);

void TG_sample_local_shrink(arma::vec& local_shrink,
                            const arma::vec& param_vec,
                            const arma::vec& prior_var_til,
                            double global_shrink,
                            double c,
                            double a);

double TG_sample_d2(double global_shrink,
                    double a,
                    double c);


#endif
