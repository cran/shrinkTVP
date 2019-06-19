#ifndef LOG_RATIO_VALUE_H
#define LOG_RATIO_VALUE_H

#include <RcppArmadillo.h>

double log_ratio_value_marginalBFS(int d, double proposal, double old_val, double scale_par, arma::vec param_vec, double b1, double b2);

#endif
