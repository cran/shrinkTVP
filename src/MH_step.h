#ifndef MH_STEP_H
#define MH_STEP_H

#include <RcppArmadillo.h>

double MH_step(double current_val, double c_tuning_par, int d, double scale_par, arma::vec param_vec, double b, double nu,
               double hyp1, double hyp2);
#endif
