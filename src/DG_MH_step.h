#ifndef DG_MH_STEP_H
#define DG_MH_STEP_H

double DG_MH_step(double current_val,
                  double a_tuning_par,
                  double scale_par,
                  const arma::vec& param_vec,
                  double b,
                  double nu,
                  bool adaptive,
                  arma::vec& batch,
                  double& curr_sd,
                  double target_rate,
                  double max_adapt,
                  int& batch_nr,
                  int batch_size,
                  int& batch_pos);

#endif
