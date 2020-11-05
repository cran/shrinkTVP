#ifndef DG_SAMPLING_FUNCTIONS_H
#define DG_SAMPLING_FUNCTIONS_H

void DG_sample_local_shrink(arma::vec& local_shrink,
                            const arma::vec& param_vec,
                            double global_shrink,
                            double a);

double DG_sample_global_shrink(const arma::vec& prior_var,
                               double a,
                               double hyper1,
                               double hyper2);

#endif
