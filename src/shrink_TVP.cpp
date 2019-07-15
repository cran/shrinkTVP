// [[Rcpp::depends(RcppArmadillo, RcppProgress)]]
#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <math.h>
#include "sample_beta_McCausland.h"
#include "sample_parameters.h"
#include "MH_step.h"
using namespace Rcpp;

// [[Rcpp::export]]
List do_shrinkTVP(arma::vec y,
                  arma::mat x,
                  arma::vec a0,
                  int niter,
                  int nburn,
                  int nthin,
                  double c0,
                  double g0,
                  double G0,
                  double d1,
                  double d2,
                  double e1,
                  double e2,
                  bool learn_lambda2,
                  bool learn_kappa2,
                  double lambda2,
                  double kappa2,
                  bool learn_a_xi,
                  bool learn_a_tau,
                  double a_xi,
                  double a_tau,
                  double c_tuning_par_xi,
                  double c_tuning_par_tau,
                  double b_xi,
                  double b_tau,
                  double nu_xi,
                  double nu_tau,
                  bool display_progress,
                  bool ret_beta_nc,
                  bool store_burn,
                  bool sv,
                  double Bsigma_sv,
                  double a0_sv,
                  double b0_sv,
                  double bmu,
                  double Bmu,
                  double y_test,
                  arma::vec x_test,
                  bool LPDS) {

  // progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // Import Rs chol function
  Environment base = Environment("package:base");
  Function Rchol = base["chol"];

  // Some necessary dimensions
  int N = y.n_elem;
  int d = x.n_cols;
  int nsave;
  if (store_burn){
    nsave = std::floor(niter/nthin);
  } else {
    nsave = std::floor((niter - nburn)/nthin);
  }

  //Storage objects
  arma::cube beta_save(N+1, d, nsave, arma::fill::none);
  arma::cube sig2_save(N,1, nsave, arma::fill::none);
  arma::mat theta_sr_save(d, nsave, arma::fill::none);
  arma::mat beta_mean_save(d, nsave, arma::fill::none);
  arma::mat xi2_save(d, nsave, arma::fill::none);
  arma::mat tau2_save(d, nsave, arma::fill::none);

  // conditional storage objects
  arma::vec kappa2_save;
  arma::vec lambda2_save;
  if (learn_kappa2){
    kappa2_save = arma::vec(nsave, arma::fill::none);
  }
  if (learn_lambda2){
    lambda2_save = arma::vec(nsave, arma::fill::none);
  }

  arma::cube beta_nc_save;
  if (ret_beta_nc){
    beta_nc_save = arma::cube(N+1, d, nsave, arma::fill::none);
  }

  arma::vec a_xi_save;
  arma::vec a_tau_save;

  int accept_a_xi_tot = 0;
  int accept_a_xi_pre = 0;
  int accept_a_xi_post = 0;

  int accept_a_tau_tot = 0;
  int accept_a_tau_pre = 0;
  int accept_a_tau_post = 0;

  if (learn_a_xi){
    a_xi_save = arma::vec(nsave, arma::fill::none);
  }
  if (learn_a_tau){
    a_tau_save = arma::vec(nsave, arma::fill::none);
  }

  arma::vec C0_save;
  arma::vec sv_mu_save;
  arma::vec sv_phi_save;
  arma::vec sv_sigma2_save;

  if (sv == false){
    C0_save = arma::vec(nsave, arma::fill::none);
  } else {
    sv_mu_save = arma::vec(nsave, arma::fill::none);
    sv_phi_save = arma::vec(nsave, arma::fill::none);
    sv_sigma2_save = arma::vec(nsave, arma::fill::none);
  }

  // Initial values and objects
  arma::mat beta_nc_samp(d, N+1, arma::fill::none);

  arma::vec beta_mean_samp(d);
  beta_mean_samp.fill(0.1);

  arma::vec theta_sr_samp(d);
  theta_sr_samp.fill(0.2);

  arma::vec tau2_samp(d);
  tau2_samp.fill(0.1);

  arma::vec xi2_samp(d);
  xi2_samp.fill(0.1);

  arma::vec xi_tau_samp = arma::join_cols(xi2_samp, tau2_samp);

  double kappa2_samp = 20;

  double lambda2_samp = 0.1;

  double a_xi_samp = 0.1;

  double a_tau_samp = 0.1;

  arma::vec kappa2_lambda_samp = {kappa2_samp, lambda2_samp};

  arma::vec h_samp(N, arma::fill::zeros);

  arma::vec alpha_samp(2*d, arma::fill::ones);

  arma::vec sig2_samp = arma::exp(h_samp);

  double C0_samp = 1;

  // SV quantities
  arma::vec sv_para = {-10, 0.5, 1};
  arma::mat mixprob(10, N);
  arma::vec mixprob_vec(mixprob.begin(), mixprob.n_elem, false);
  arma::ivec r(N);
  double h0 = -10;
  double B011inv         = 1e-8;
  double B022inv         = 1e-12;
  bool Gammaprior        = true;
  double MHcontrol       = -1;
  int parameterization   = 3;
  bool centered_baseline = parameterization % 2; // 1 for C, 0 for NC baseline
  int MHsteps = 2;
  bool dontupdatemu = 0;
  double cT = N/2.0;
  double C0_sv = 1.5*Bsigma_sv;
  bool truncnormal = false;
  double priorlatent0 = -1;

  // Values for LPDS
  arma::vec m_N_samp;
  arma::mat chol_C_N_inv_samp;
  double sum_pred_y = 0;
  arma::vec F;
  arma::vec LF;
  double S;
  double y_hat;
  double sig2_pred;

  // Override inital values with user specified fixed values
  if (learn_kappa2 == false){
    kappa2_samp = kappa2;
  }
  if (learn_lambda2 == false){
    lambda2_samp = lambda2;
  }

  if (!learn_a_xi){
    a_xi_samp = a_xi;
  }
  if (!learn_a_tau){
    a_tau_samp = a_tau;
  }

  // Values to check if the sampler failed or not
  bool succesful = true;
  std::string fail;
  int fail_iter;


  // Introduce additional index post_j that is used to calculate accurate storage positions in case of thinning
  int post_j = 1;

  // Begin Gibbs loop
  for (int j = 0; j < niter; j++){
    // step a)
    // sample time varying beta.tilde parameters (NC parametrization)
    try {
      sample_beta_McCausland(beta_nc_samp, y, x, theta_sr_samp, sig2_samp, beta_mean_samp, m_N_samp, chol_C_N_inv_samp, LPDS, N, d, Rchol);
    } catch (...){
      beta_nc_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample beta_nc";
        fail_iter = j + 1;
        succesful = false;
      }
    }


    // step b)
    // sample alpha
    arma::mat x_tilde = x % (beta_nc_samp.cols(1,N)).t();
    arma::mat W = arma::join_rows(x, x_tilde);

    try {
      sample_alpha(alpha_samp, y, x, W, tau2_samp, xi2_samp, sig2_samp, a0, d, Rchol);
    } catch(...){
      alpha_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample alpha";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // Weave back into centered parameterization
    beta_mean_samp = alpha_samp.rows(0, d-1);
    theta_sr_samp = alpha_samp.rows(d, 2*d-1);
    arma::mat beta_nc_samp_tilde = beta_nc_samp.each_col() % theta_sr_samp;
    arma::mat betaenter = beta_nc_samp_tilde.each_col() + beta_mean_samp;

    // Difference beta outside of function (for numerical stability)
    arma::mat beta_diff_pre = arma::diff(beta_nc_samp, 1, 1);
    arma::mat beta_diff =  beta_diff_pre.each_col() % theta_sr_samp;

    // step c)
    // resample alpha
    try {
      resample_alpha_diff(alpha_samp, betaenter, theta_sr_samp, beta_mean_samp, beta_diff, xi2_samp, tau2_samp, d, N);
    } catch(...) {
      alpha_samp.fill(nanl(""));
      if (succesful == true){
        fail = "resample alpha";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // Calculate NC betas with new alpha
    beta_mean_samp = alpha_samp.rows(0, d-1);
    theta_sr_samp = alpha_samp.rows(d, 2*d-1);
    beta_nc_samp = betaenter.each_col() - beta_mean_samp;
    beta_nc_samp.each_col() /= theta_sr_samp;

    x_tilde = x % (beta_nc_samp.cols(1,N)).t();
    W = arma::join_rows(x, x_tilde);

    // step d)
    // sample a_xi and a_tau with MH
    if (learn_a_xi){
      double before = a_xi_samp;
      try {
        a_xi_samp = MH_step(a_xi_samp, c_tuning_par_xi, d, kappa2_samp, theta_sr_samp, b_xi , nu_xi, d1, d2);
      } catch(...) {
        a_xi_samp = nanl("");
        if (succesful == true){
          fail = "sample a_xi";
          fail_iter = j + 1;
          succesful = false;
        }
      }

      if (before != a_xi_samp){
        accept_a_xi_tot += 1;
        if (j < nburn){
          accept_a_xi_pre += 1;
        } else {
          accept_a_xi_post += 1;
        }
      }
    }


    if (learn_a_tau){
      double before = a_tau_samp;
      try {
        a_tau_samp = MH_step(a_tau_samp, c_tuning_par_tau, d, lambda2_samp, beta_mean_samp, b_tau , nu_tau, e1, e2);
      } catch(...){
        a_tau_samp = nanl("");
        if (succesful == true){
          fail = "sample a_tau";
          fail_iter = j + 1;
          succesful = false;
        }
      }

      if (before != a_tau_samp){
        accept_a_tau_tot += 1;
        if (j < nburn){
          accept_a_tau_pre += 1;
        } else {
          accept_a_tau_post += 1;
        }
      }
    }

    // step e)
    // sample tau2 and xi2
    try {
      sample_tau2(tau2_samp, beta_mean_samp, lambda2_samp, a_tau_samp, d);
    } catch(...) {
      tau2_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample tau2";
        fail_iter = j + 1;
        succesful = false;
      }
    }
    try {
      sample_xi2(xi2_samp, theta_sr_samp, kappa2_samp, a_xi_samp, d);
    } catch(...) {
      xi2_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample xi2";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // sample kappa2 and lambda2, if the user specified it
    if (learn_kappa2){
      try {
        kappa2_samp = sample_kappa2(xi2_samp, a_xi_samp, d1, d2, d);
      } catch (...) {
        kappa2_samp = nanl("");
        if (succesful == true){
          fail = "sample kappa2";
          fail_iter = j + 1;
          succesful = false;
        }
      }
    }

    if (learn_lambda2){
      try {
        lambda2_samp = sample_lambda2(tau2_samp, a_tau_samp, e1, e2, d);
      } catch (...) {
        lambda2_samp = nanl("");
        if (succesful == true){
          fail = "sample lambda2";
          fail_iter = j + 1;
          succesful = false;
        }
      }
    }

    // step f)
    // sample sigma2 from homoscedastic or SV case
    try {
      if (sv){
        arma::vec datastand = arma::log(arma::square(y - x * beta_mean_samp - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp));

        arma::vec cur_h = arma::log(sig2_samp);
        stochvol::update_sv(datastand, sv_para, cur_h, h0, mixprob_vec, r, centered_baseline, C0_sv, cT,
                            Bsigma_sv, a0_sv, b0_sv, bmu, Bmu, B011inv, B022inv, Gammaprior,
                            truncnormal, MHcontrol, MHsteps, parameterization, dontupdatemu, priorlatent0);

        sig2_samp = arma::exp(cur_h);
      } else {
        sample_sigma2(sig2_samp, y, W, alpha_samp, c0, C0_samp, N);
      }
    } catch(...) {
      sig2_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample sigma2";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    if(sv == false){
      try {
        C0_samp = sample_C0(sig2_samp, g0, c0, G0);
      } catch(...) {
        C0_samp = nanl("");
        if (succesful == true){
          fail = "sample C0";
          fail_iter = j + 1;
          succesful = false;
        }
      }
    }

    // adjust start of storage point depending on store_burn
    int nburn_new = nburn;
    if(store_burn){
      nburn_new = 0;
    }

    // Calculate predicted value and add to sum for LPDS
    if((LPDS == true) && (j % nthin == 0) && (j >= nburn)){

      if (sv == false){
        sig2_pred = sig2_samp(0);
      } else {
        sig2_pred = std::exp(R::rnorm(sv_para(0) + sv_para(1) * (std::log(sig2_samp(N-1)) - sv_para(0)), std::sqrt(sv_para(2))));
      }

      F = x_test % theta_sr_samp;
      LF = arma::solve(arma::trimatl(chol_C_N_inv_samp), F);
      S = arma::as_scalar(LF.t() * LF + F.t() * F + sig2_pred);
      y_hat = arma::as_scalar(x_test.t() * beta_mean_samp + F.t() * m_N_samp);
      sum_pred_y += R::dnorm(y_test, y_hat, std::sqrt(S), 0);
    }


    // Increment index i if burn-in period is over
    if (j > nburn_new){
      post_j++;
    }

    // Store everything
    if ((post_j % nthin == 0) && (j >= nburn_new)){
      // Caluclate beta
      // This is in the if condition to save unnecessary computations if beta is not saved
      arma::mat beta =  (beta_nc_samp.each_col() % theta_sr_samp).each_col() + beta_mean_samp;

      sig2_save.slice((post_j-1)/nthin) = sig2_samp;
      theta_sr_save.col((post_j-1)/nthin) = theta_sr_samp;
      beta_mean_save.col((post_j-1)/nthin) = beta_mean_samp;
      beta_save.slice((post_j-1)/nthin) = beta.t();
      xi2_save.col((post_j-1)/nthin) = xi2_samp;
      tau2_save.col((post_j-1)/nthin) = tau2_samp;

      //conditional storing
      if (ret_beta_nc){
        beta_nc_save.slice((post_j-1)/nthin) = beta_nc_samp.t();
      }

      if (learn_kappa2){
        kappa2_save((post_j-1)/nthin) = kappa2_samp;
      }
      if (learn_lambda2){
        lambda2_save((post_j-1)/nthin) = lambda2_samp;
      }

      if (learn_a_xi){
        a_xi_save((post_j-1)/nthin) = a_xi_samp;
      }

      if (learn_a_tau){
        a_tau_save((post_j-1)/nthin) = a_tau_samp;
      }

      if (sv == false){
        C0_save((post_j-1)/nthin) = C0_samp;
      } else {
        sv_mu_save((post_j-1)/nthin) = sv_para(0);
        sv_phi_save((post_j-1)/nthin) = sv_para(1);
        sv_sigma2_save((post_j-1)/nthin) = sv_para(2);
      }
    }

    // Random sign switch
    for (int i = 0; i < d; i++){
      if(R::runif(0,1) > 0.5){
        theta_sr_samp(i) = -theta_sr_samp(i);
      }
    }

    // Increment progress bar
    if (arma::any(prog_rep_points == j)){
      p.increment();
    }

    // Check for user interrupts
    if (j % 500 == 0){
      Rcpp::checkUserInterrupt();
    }

    // Break loop if succesful is false
    if (!succesful){
      break;
    }
  }

  arma::mat LPDS_res(1,1);
  LPDS_res(0,0) = std::log(sum_pred_y/std::floor((niter - nburn)/nthin));

  // return everything as a nested list (due to size restrictions on Rcpp::Lists)
  return List::create(_["sigma2"] = sig2_save,
                      _["theta_sr"] = theta_sr_save.t(),
                      _["beta_mean"] = beta_mean_save.t(),
                      _["beta_nc"] = beta_nc_save,
                      _["beta"] = beta_save,
                      _["xi2"] = xi2_save.t(),
                      _["a_xi"] = a_xi_save,
                      _["a_xi_acceptance"] = List::create(
                        _["a_xi_acceptance_total"] = (double)accept_a_xi_tot/niter,
                        _["a_xi_acceptance_pre"] = (double)accept_a_xi_pre/nburn,
                        _["a_xi_acceptance_post"] = (double)accept_a_xi_post/(niter - nburn)),
                      _["tau2"] = tau2_save.t(),
                      _["a_tau"] = a_tau_save,
                      _["a_tau_acceptance"] = List::create(
                          _["a_tau_acceptance_total"] = (double)accept_a_tau_tot/niter,
                          _["a_tau_acceptance_pre"] = (double)accept_a_tau_pre/nburn,
                          _["a_tau_acceptance_post"] = (double)accept_a_tau_post/(niter - nburn)),
                      _["kappa2"] = kappa2_save,
                      _["lambda2"] = lambda2_save,
                      _["C0"] = C0_save,
                      _["sv_mu"] = sv_mu_save,
                      _["sv_phi"] = sv_phi_save,
                      _["sv_sigma2"] = sv_sigma2_save,
                      _["LPDS"] = LPDS_res,
                      _["success_vals"] = List::create(
                        _["success"] = succesful,
                        _["fail"] = fail,
                        _["fail_iter"] = fail_iter
                      ));
}



