#ifndef RHO_P_MH_STEP_H
#define RHO_P_MH_STEP_H


double rho_p_MH_step_marg_oeverything(double current_val,
                                      const arma::vec& psi,
                                      double psi_0,
                                      double a,
                                      double c,
                                      double a_rho,
                                      double b_rho,
                                      double alpha,
                                      double beta,
                                      double rho_tuning_par,
                                      bool adaptive,
                                      arma::vec& batch,
                                      double& curr_sd,
                                      double target_rate,
                                      double max_adapt,
                                      int& batch_nr,
                                      int batch_size,
                                      int& batch_pos);

#endif

