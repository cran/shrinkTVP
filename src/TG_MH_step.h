#ifndef TG_MH_STEP_H
#define TG_MH_STEP_H

double TG_MH_step(double current_val,
                  double c_tuning_par,
                  double scale_par,
                  const arma::vec& scale_vec,
                  const arma::vec& param_vec,
                  double b,
                  double nu,
                  bool is_c,
                  double scale_scale,
                  double other_hyp,
                  bool common_shrink_par,
                  bool adaptive,
                  arma::vec& batch,
                  double& curr_sd,
                  double target_rate,
                  double max_adapt,
                  int& batch_nr,
                  int batch_size,
                  int& batch_pos);

#endif
