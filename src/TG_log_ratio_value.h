#ifndef TG_LOG_RATIO_VALUE_H
#define TG_LOG_RATIO_VALUE_H

double TG_log_ratio_value_marginalBFS(double proposal,
                                      double old_val,
                                      double scale_par,
                                      const arma::vec& scale_vec,
                                      const arma::vec& param_vec,
                                      double scale_scale,
                                      double c,
                                      double b1,
                                      double b2,
                                      bool common_shrink_par);

double TG_log_ratio_value_tg(double proposal,
                             double old_val,
                             double scale_par,
                             const arma::vec& scale_vec,
                             const arma::vec& param_vec,
                             double scale_scale,
                             double a,
                             double b1,
                             double b2);

#endif



