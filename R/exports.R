pred_dens_mix_approx_dyn <- function(x_test, y_test, theta_sr, beta_mean, sig2_samp, psi_future, sv, sv_phi, sv_mu, sv_sigma2, chol_C_N_inv_samp, m_N_samp, M, log) {
  .Call(`_shrinkTVP_pred_dens_mix_approx_dyn`, x_test, y_test, theta_sr, beta_mean, sig2_samp, psi_future, sv, sv_phi, sv_mu, sv_sigma2, chol_C_N_inv_samp, m_N_samp, M, log)
}

calc_fitted_cpp <- function(y, x, beta) {
    .Call(`_shrinkTVP_calc_fitted_cpp`, y, x, beta)
}

shrinkTVP_cpp <- function(y, x, mod_type, niter, nburn, nthin, c0, g0, G0, d1, d2, e1, e2, learn_lambda2_B, learn_kappa2_B, lambda2_B, kappa2_B, learn_a_xi, learn_a_tau, a_xi, a_tau, learn_c_xi, learn_c_tau, c_xi, c_tau, a_eq_c_xi, a_eq_c_tau, a_tuning_par_xi, a_tuning_par_tau, c_tuning_par_xi, c_tuning_par_tau, beta_a_xi, beta_a_tau, alpha_a_xi, alpha_a_tau, beta_c_xi, beta_c_tau, alpha_c_xi, alpha_c_tau, display_progress, sv, Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, adaptive, target_rates, max_adapts, batch_sizes, starting_vals) {
    .Call(`_shrinkTVP_shrinkTVP_cpp`, y, x, mod_type, niter, nburn, nthin, c0, g0, G0, d1, d2, e1, e2, learn_lambda2_B, learn_kappa2_B, lambda2_B, kappa2_B, learn_a_xi, learn_a_tau, a_xi, a_tau, learn_c_xi, learn_c_tau, c_xi, c_tau, a_eq_c_xi, a_eq_c_tau, a_tuning_par_xi, a_tuning_par_tau, c_tuning_par_xi, c_tuning_par_tau, beta_a_xi, beta_a_tau, alpha_a_xi, alpha_a_tau, beta_c_xi, beta_c_tau, alpha_c_xi, alpha_c_tau, display_progress, sv, Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, adaptive, target_rates, max_adapts, batch_sizes, starting_vals)
}

shrinkDTVP_cpp <- function(y, x, mod_type, iid, niter, nburn, nthin, c0, g0, G0, d1, d2, e1, e2, learn_lambda2_B, learn_kappa2_B, lambda2_B, kappa2_B, learn_a_xi, learn_a_tau, a_xi, a_tau, learn_c_xi, learn_c_tau, c_xi, c_tau, a_eq_c_xi, a_eq_c_tau, a_tuning_par_xi, a_tuning_par_tau, c_tuning_par_xi, c_tuning_par_tau, beta_a_xi, beta_a_tau, alpha_a_xi, alpha_a_tau, beta_c_xi, beta_c_tau, alpha_c_xi, alpha_c_tau, alpha_rho, beta_rho, a_psi, c_psi, a_rho, b_rho, shrink_inter, display_progress, sv, Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, adaptive_rho, tuning_par_rho, target_rate_rho, max_adapt_rho, batch_size_rho, adaptive, target_rates, max_adapts, batch_sizes, starting_vals) {
  .Call(`_shrinkTVP_shrinkDTVP_cpp`, y, x, mod_type, iid, niter, nburn, nthin, c0, g0, G0, d1, d2, e1, e2, learn_lambda2_B, learn_kappa2_B, lambda2_B, kappa2_B, learn_a_xi, learn_a_tau, a_xi, a_tau, learn_c_xi, learn_c_tau, c_xi, c_tau, a_eq_c_xi, a_eq_c_tau, a_tuning_par_xi, a_tuning_par_tau, c_tuning_par_xi, c_tuning_par_tau, beta_a_xi, beta_a_tau, alpha_a_xi, alpha_a_tau, beta_c_xi, beta_c_tau, alpha_c_xi, alpha_c_tau, alpha_rho, beta_rho, a_psi, c_psi, a_rho, b_rho, shrink_inter, display_progress, sv, Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, adaptive_rho, tuning_par_rho, target_rate_rho, max_adapt_rho, batch_size_rho, adaptive, target_rates, max_adapts, batch_sizes, starting_vals)
}




methods::setLoadAction(function(ns) {
  .Call('_shrinkTVP_RcppExport_registerCCallable', PACKAGE = 'shrinkTVP')
})

