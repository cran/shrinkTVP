#ifndef DG_LOG_RATIO_VALUE_H
#define DG_LOG_RATIO_VALUE_H

double DG_log_ratio_value_marginalBFS(double proposal,
                                      double old_val,
                                      double scale_par,
                                      const arma::vec& param_vec,
                                      double b1,
                                      double b2);

#endif



