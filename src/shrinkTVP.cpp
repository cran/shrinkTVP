#include <RcppArmadillo.h>
#include <stochvol.h>
#include <progress.hpp>
#include <math.h>
#include "sample_beta_McCausland.h"
#include "common_sampling_functions.h"
#include "sample_TG_TVP.h"
#include "sample_DG_TVP.h"
#include "cpp_utilities.h"
using namespace Rcpp;



List shrinkTVP_cpp(arma::vec y,
                   arma::mat x,
                   std::string mod_type,
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
                   bool learn_lambda2_B,
                   bool learn_kappa2_B,
                   double lambda2_B,
                   double kappa2_B,
                   bool learn_a_xi,
                   bool learn_a_tau,
                   double a_xi,
                   double a_tau,
                   bool learn_c_xi,
                   bool learn_c_tau,
                   double c_xi,
                   double c_tau,
                   bool a_eq_c_xi,
                   bool a_eq_c_tau,
                   double a_tuning_par_xi,
                   double a_tuning_par_tau,
                   double c_tuning_par_xi,
                   double c_tuning_par_tau,
                   double beta_a_xi,
                   double beta_a_tau,
                   double alpha_a_xi,
                   double alpha_a_tau,
                   double beta_c_xi,
                   double beta_c_tau,
                   double alpha_c_xi,
                   double alpha_c_tau,
                   bool display_progress,
                   bool sv,
                   double Bsigma_sv,
                   double a0_sv,
                   double b0_sv,
                   double bmu,
                   double Bmu,
                   arma::vec adaptive,
                   arma::vec target_rates,
                   arma::vec max_adapts,
                   arma::ivec batch_sizes,
                   Rcpp::List starting_vals) {

  // Progress bar setup
  arma::vec prog_rep_points = arma::round(arma::linspace(0, niter, 50));
  Progress p(50, display_progress);

  // Define some necessary dimensions
  int N = y.n_elem;
  int d = x.n_cols;
  int nsave = std::floor((niter - nburn)/nthin);

  // Storage objects (returned to R at the end)
  arma::cube beta_save(N+1, d, nsave, arma::fill::zeros);
  arma::cube sigma2_save(N,1, nsave, arma::fill::zeros);
  arma::mat theta_sr_save(d, nsave, arma::fill::zeros);
  arma::mat beta_mean_save(d, nsave, arma::fill::zeros);

  // conditional storage objects (only used for some model types)
  arma::mat xi2_save;
  arma::mat tau2_save;
  if (mod_type != "ridge") {
    xi2_save = arma::mat(d, nsave, arma::fill::zeros);
    tau2_save = arma::mat(d, nsave, arma::fill::zeros);
  }

  arma::mat lambda2_save;
  arma::mat kappa2_save;
  if (mod_type == "triple") {
    lambda2_save = arma::mat(d, nsave, arma::fill::zeros);
    kappa2_save = arma::mat(d, nsave, arma::fill::zeros);
  }

  arma::vec kappa2_B_save;
  arma::vec lambda2_B_save;
  if (learn_kappa2_B & (mod_type != "ridge")) {
    kappa2_B_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_lambda2_B & (mod_type != "ridge")) {
    lambda2_B_save = arma::vec(nsave, arma::fill::zeros);
  }

  arma::vec a_xi_save;
  arma::vec a_tau_save;
  arma::vec c_xi_save;
  arma::vec c_tau_save;

  if (learn_a_xi & (mod_type != "ridge")) {
    a_xi_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_a_tau & (mod_type != "ridge")) {
    a_tau_save = arma::vec(nsave, arma::fill::zeros);
  }

  if (learn_c_xi & (mod_type == "triple")) {
    c_xi_save = arma::vec(nsave, arma::fill::zeros);
  }
  if (learn_c_tau & (mod_type == "triple")) {
    c_tau_save = arma::vec(nsave, arma::fill::zeros);
  }

  arma::vec C0_save;
  arma::vec sv_mu_save;
  arma::vec sv_phi_save;
  arma::vec sv_sigma2_save;

  if (sv == false) {
    C0_save = arma::vec(nsave, arma::fill::zeros);
  } else {
    sv_mu_save = arma::vec(nsave, arma::fill::zeros);
    sv_phi_save = arma::vec(nsave, arma::fill::zeros);
    sv_sigma2_save = arma::vec(nsave, arma::fill::zeros);
  }

  arma::vec a_xi_sd_save;
  arma::vec a_tau_sd_save;
  arma::vec c_xi_sd_save;
  arma::vec c_tau_sd_save;
  arma::vec a_xi_acc_rate_save;
  arma::vec a_tau_acc_rate_save;
  arma::vec c_xi_acc_rate_save;
  arma::vec c_tau_acc_rate_save;

  if (bool(adaptive(0)) & learn_a_xi & (mod_type != "ridge")) {
    a_xi_sd_save = arma::vec(std::floor(niter/batch_sizes(0)), arma::fill::zeros);
    a_xi_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(0)), arma::fill::zeros);
  }

  if (bool(adaptive(1)) & learn_a_tau & (mod_type != "ridge")) {
    a_tau_sd_save = arma::vec(std::floor(niter/batch_sizes(1)), arma::fill::zeros);
    a_tau_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(1)), arma::fill::zeros);
  }

  if (bool(adaptive(2)) & learn_c_xi & (mod_type == "triple")) {
    c_xi_sd_save = arma::vec(std::floor(niter/batch_sizes(2)), arma::fill::zeros);
    c_xi_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(2)), arma::fill::zeros);
  }

  if (bool(adaptive(3)) & learn_c_tau & (mod_type == "triple")) {
    c_tau_sd_save = arma::vec(std::floor(niter/batch_sizes(3)), arma::fill::zeros);
    c_tau_acc_rate_save = arma::vec(std::floor(niter/batch_sizes(3)), arma::fill::zeros);
  }

  // Initial values and objects holding current values of samples
  // Parameters where learning can be toggled are also overwritten if learning is not activated
  arma::mat beta_samp(d, N+1, arma::fill::zeros);
  arma::mat beta_nc_samp(d, N+1, arma::fill::zeros);

  arma::vec beta_mean_samp(d);
  beta_mean_samp = as<arma::vec>(starting_vals["beta_mean_st"]);

  arma::vec theta_sr_samp(d);
  theta_sr_samp = as<arma::vec>(starting_vals["theta_sr_st"]);

  // Conditional objects
  double kappa2_B_samp;
  if (learn_kappa2_B & (mod_type != "ridge")) {
    kappa2_B_samp = as<double>(starting_vals["kappa2_B_st"]);
  } else {
    kappa2_B_samp = kappa2_B;
  }

  double lambda2_B_samp;
  if (learn_lambda2_B & (mod_type != "ridge")) {
    lambda2_B_samp = as<double>(starting_vals["lambda2_B_st"]);
  } else {
    lambda2_B_samp = lambda2_B;
  }

  double a_xi_samp;
  if (learn_a_xi & (mod_type != "ridge")) {
    a_xi_samp = as<double>(starting_vals["a_xi_st"]);
  } else {
    a_xi_samp = a_xi;
  }

  double a_tau_samp;
  if (learn_a_tau & (mod_type != "ridge")) {
    a_tau_samp = as<double>(starting_vals["a_tau_st"]);
  } else {
    a_tau_samp = a_tau;
  }

  double c_xi_samp;
  if (learn_c_xi & (mod_type == "triple")) {
    c_xi_samp = as<double>(starting_vals["c_xi_st"]);
  } else {
    c_xi_samp = c_xi;
  }

  double c_tau_samp;
  if (learn_c_tau & (mod_type == "triple")) {
    c_tau_samp = as<double>(starting_vals["c_tau_st"]);
  } else {
    c_tau_samp = c_tau;
  }

  double d2_samp;
  double e2_samp;

  arma::vec lambda2_til_samp(d);
  arma::vec kappa2_til_samp(d);
  if (mod_type == "triple") {
    lambda2_til_samp = as<arma::vec>(starting_vals["lambda2_st"]);
    kappa2_til_samp = as<arma::vec>(starting_vals["kappa2_st"]);

  }

  arma::vec tau2_samp(d);
  arma::vec xi2_samp(d);
  arma::vec tau2_til_samp(d);
  arma::vec xi2_til_samp(d);
  if (mod_type == "double") {
    tau2_samp = as<arma::vec>(starting_vals["tau2_st"]);
    xi2_samp = as<arma::vec>(starting_vals["xi2_st"]);
  } else if (mod_type == "triple") {
    tau2_til_samp = as<arma::vec>(starting_vals["tau2_st"]);
    xi2_til_samp = as<arma::vec>(starting_vals["xi2_st"]);

    calc_xi2_tau2(xi2_samp,
                  xi2_til_samp,
                  kappa2_til_samp,
                  kappa2_B_samp,
                  c_xi_samp,
                  a_xi_samp);

    calc_xi2_tau2(tau2_samp,
                  tau2_til_samp,
                  lambda2_til_samp,
                  lambda2_B_samp,
                  c_tau_samp,
                  a_tau_samp);
  } else {
    tau2_samp.fill(2.0/lambda2_B_samp);
    xi2_samp.fill(2.0/kappa2_B_samp);
  }


  arma::vec sigma2_samp;
  double C0_samp;
  arma::vec sv_para(3);
  if (sv == true) {
    sigma2_samp = as<arma::vec>(starting_vals["sigma2_st"]);
    sv_para = {as<double>(starting_vals["sv_mu_st"]),
               as<double>(starting_vals["sv_phi_st"]),
               as<double>(starting_vals["sv_sigma2_st"])};
  } else {
    sigma2_samp.copy_size(y);
    sigma2_samp.fill(as<double>(starting_vals["sigma2_st"]));
    C0_samp = as<double>(starting_vals["C0_st"]);
  }

  // Objects required for stochvol to work
  arma::uvec r(N); r.fill(5);
  double h0_samp = as<double>(starting_vals["h0_st"]);
  using stochvol::PriorSpec;
  const PriorSpec prior_spec = {  // prior specification object for the update_*_sv functions
    PriorSpec::Latent0(),  // stationary prior distribution on priorlatent0
    PriorSpec::Mu(PriorSpec::Normal(bmu, std::sqrt(Bmu))),  // normal prior on mu
    PriorSpec::Phi(PriorSpec::Beta(a0_sv, b0_sv)),  // stretched beta prior on phi
    PriorSpec::Sigma2(PriorSpec::Gamma(0.5, 0.5 / Bsigma_sv))  // normal(0, Bsigma) prior on sigma
  };  // heavy-tailed, leverage, regression turned off
  using stochvol::ExpertSpec_FastSV;
  const ExpertSpec_FastSV expert {  // very expert settings for the Kastner, Fruehwirth-Schnatter (2014) sampler
    true,  // interweave
    stochvol::Parameterization::CENTERED,  // centered baseline always
    1e-8,  // B011inv,
    1e-12,  //B022inv,
    2,  // MHsteps,
    ExpertSpec_FastSV::ProposalSigma2::INDEPENDENCE,  // independece proposal for sigma
    -1,  // unused for independence prior for sigma
    ExpertSpec_FastSV::ProposalPhi::IMMEDIATE_ACCEPT_REJECT_NORMAL  // immediately reject (mu,phi,sigma) if proposed phi is outside (-1, 1)
  };

  // Objects necessary for the use of adaptive MH
  arma::mat batches;
  arma::vec curr_sds;
  arma::ivec batch_nrs;
  arma::ivec batch_pos;

  // This ensures that adaptive MH is run correctly in the standalone case
  try {
    Rcpp::List temp = starting_vals["internals"];
    batch_pos = as<arma::ivec>(temp["batch_pos_st"]);
    batches = as<arma::mat>(temp["batches_st"]);
    curr_sds = as<arma::vec>(temp["curr_sds_st"]);
    batch_nrs = as<arma::ivec>(temp["batch_nrs_st"]);
  } catch(...) {
    batch_pos = arma::ivec(4, arma::fill::zeros);
    batches = arma::mat(arma::max(batch_sizes), 4, arma::fill::zeros);
    curr_sds = {a_tuning_par_xi,
                a_tuning_par_tau,
                c_tuning_par_xi,
                c_tuning_par_tau};
    batch_nrs = arma::ivec(4, arma::fill::ones);
  }

  // Values for LPDS calculation
  arma::cube m_N_save(d, 1, nsave);
  arma::cube chol_C_N_inv_save(d, d, nsave);
  arma::vec m_N_samp;
  arma::mat chol_C_N_inv_samp;

  // Values to check if the sampler failed or not
  bool succesful = true;
  std::string fail;
  int fail_iter;


  // Introduce additional index post_j that is used to calculate accurate storage positions in case of thinning
  int post_j = 1;
  // Begin Gibbs loop
  for (int j = 0; j < niter; j++) {
    // sample time varying beta.tilde parameters (NC parametrization)
    try {
      sample_beta_McCausland(beta_nc_samp,
                             y,
                             x,
                             theta_sr_samp,
                             sigma2_samp,
                             beta_mean_samp,
                             m_N_samp,
                             chol_C_N_inv_samp);
    } catch (...) {
      beta_nc_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "sample beta_nc";
        fail_iter = j + 1;
        succesful = false;
      }
    }


    // Sample alpha (theta_sr and beta_mean)
    try {
      sample_alpha(beta_mean_samp,
                   theta_sr_samp,
                   y,
                   x,
                   beta_nc_samp,
                   sigma2_samp,
                   tau2_samp,
                   xi2_samp);
    } catch(...){
      beta_mean_samp.fill(nanl(""));
      theta_sr_samp.fill(nanl(""));
      if (succesful == true){
        fail = "sample alpha";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // Weave into centered parameterization and resample alpha
    to_CP(beta_samp,
          beta_nc_samp,
          theta_sr_samp,
          beta_mean_samp);

    try {
      resample_alpha(beta_mean_samp,
                     theta_sr_samp,
                     beta_samp,
                     beta_nc_samp,
                     xi2_samp,
                     tau2_samp);
    } catch(...) {
      beta_mean_samp.fill(nanl(""));
      theta_sr_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "resample alpha";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    // Weave back into non-centered parameterization
    to_NCP(beta_nc_samp,
           beta_samp,
           theta_sr_samp,
           beta_mean_samp);

    // Sample prior variance, differentiating between ridge, double gamma and triple gamma
    if (mod_type == "double") {
      sample_DG_TVP(beta_mean_samp,
                    theta_sr_samp,
                    tau2_samp,
                    xi2_samp,
                    lambda2_B_samp,
                    kappa2_B_samp,
                    a_xi_samp,
                    beta_a_xi,
                    alpha_a_xi,
                    a_tau_samp,
                    beta_a_tau,
                    alpha_a_tau,
                    d1,
                    d2,
                    e1,
                    e2,
                    learn_kappa2_B,
                    learn_lambda2_B,
                    learn_a_xi,
                    learn_a_tau,
                    a_tuning_par_xi,
                    a_tuning_par_tau,
                    adaptive,
                    batches,
                    curr_sds,
                    target_rates,
                    max_adapts,
                    batch_nrs,
                    batch_sizes,
                    batch_pos,
                    j,
                    succesful,
                    fail,
                    fail_iter);
    } else if (mod_type == "triple") {
      sample_TG_TVP(beta_mean_samp,
                    theta_sr_samp,
                    tau2_samp,
                    xi2_samp,
                    tau2_til_samp,
                    xi2_til_samp,
                    lambda2_til_samp,
                    kappa2_til_samp,
                    lambda2_B_samp,
                    kappa2_B_samp,
                    a_xi_samp,
                    beta_a_xi,
                    alpha_a_xi,
                    a_tau_samp,
                    beta_a_tau,
                    alpha_a_tau,
                    d2_samp,
                    e2_samp,
                    c_xi_samp,
                    c_tau_samp,
                    beta_c_xi,
                    alpha_c_xi,
                    beta_c_tau,
                    alpha_c_tau,
                    learn_kappa2_B,
                    learn_lambda2_B,
                    learn_a_xi,
                    learn_a_tau,
                    learn_c_xi,
                    learn_c_tau,
                    a_tuning_par_xi,
                    a_tuning_par_tau,
                    c_tuning_par_xi,
                    c_tuning_par_tau,
                    a_eq_c_xi,
                    a_eq_c_tau,
                    adaptive,
                    batches,
                    curr_sds,
                    target_rates,
                    max_adapts,
                    batch_nrs,
                    batch_sizes,
                    batch_pos,
                    j,
                    succesful,
                    fail,
                    fail_iter);
    }

    // sample sigma2 from homoscedastic or SV case
    try {
      if (sv) {

        arma::vec datastand = 2 * arma::log(arma::abs(y - x * beta_mean_samp - (x % beta_nc_samp.cols(1,N).t()) * theta_sr_samp));
        std::for_each(datastand.begin(), datastand.end(), res_protector);

        // update_sv needs sigma and not sigma^2
        double mu = sv_para(0);
        double phi = sv_para(1);
        double sigma = std::sqrt(sv_para(2));

        arma::vec cur_h = arma::log(sigma2_samp);
        stochvol::update_fast_sv(datastand, mu, phi, sigma, h0_samp, cur_h, r, prior_spec, expert);

        // Write back into sample object
        sigma2_samp = arma::exp(cur_h);

        // change back to sigma^2
        sv_para = {mu, phi, std::pow(sigma, 2)};

        std::for_each(sigma2_samp.begin(), sigma2_samp.end(), res_protector);

      } else {
        sample_sigma2(sigma2_samp,
                      y,
                      x,
                      beta_nc_samp,
                      beta_mean_samp,
                      theta_sr_samp,
                      c0,
                      C0_samp);
      }
    } catch(...) {
      sigma2_samp.fill(nanl(""));
      if (succesful == true) {
        fail = "sample sigma2";
        fail_iter = j + 1;
        succesful = false;
      }
    }

    if(sv == false) {
      try {
        C0_samp =  sample_C0(sigma2_samp,
                             g0,
                             c0,
                             G0);
      } catch(...) {
        C0_samp = nanl("");
        if (succesful == true) {
          fail = "sample C0";
          fail_iter = j + 1;
          succesful = false;
        }
      }
    }

    // Increment index post_j if burn-in period is over
    if (j > nburn) {
      post_j++;
    }

    // Store everything (skip the steps when the sampler is used in the standalone case)
    if (!(niter == 1 && nburn == 0)){
      if ((post_j % nthin == 0) && (j >= nburn)) {
        // Caluclate beta
        // This is in the if condition to save unnecessary computations if beta is not saved
        to_CP(beta_samp,
              beta_nc_samp,
              theta_sr_samp,
              beta_mean_samp);

        sigma2_save.slice((post_j-1)/nthin) = sigma2_samp;
        theta_sr_save.col((post_j-1)/nthin) = theta_sr_samp;
        beta_mean_save.col((post_j-1)/nthin) = beta_mean_samp;
        beta_save.slice((post_j-1)/nthin) = beta_samp.t();
        m_N_save.slice((post_j-1)/nthin) = m_N_samp;
        chol_C_N_inv_save.slice((post_j-1)/nthin) = chol_C_N_inv_samp;

        // Conditional storing
        if (mod_type == "double") {
          xi2_save.col((post_j-1)/nthin) = xi2_samp;
          tau2_save.col((post_j-1)/nthin) = tau2_samp;
        } else if (mod_type == "triple") {
          xi2_save.col((post_j-1)/nthin) = xi2_til_samp;
          tau2_save.col((post_j-1)/nthin) = tau2_til_samp;
        }

        if (mod_type == "triple") {
          kappa2_save.col((post_j-1)/nthin) = kappa2_til_samp;
          lambda2_save.col((post_j-1)/nthin) = lambda2_til_samp;
        }

        if (learn_kappa2_B & (mod_type != "ridge")) {
          kappa2_B_save((post_j-1)/nthin) = kappa2_B_samp;
        }
        if (learn_lambda2_B & (mod_type != "ridge")) {
          lambda2_B_save((post_j-1)/nthin) = lambda2_B_samp;
        }

        if (learn_a_xi & (mod_type != "ridge")) {
          a_xi_save((post_j-1)/nthin) = a_xi_samp;
        }

        if (learn_a_tau & (mod_type != "ridge")) {
          a_tau_save((post_j-1)/nthin) = a_tau_samp;
        }

        if (learn_c_xi & (mod_type == "triple")) {
          c_xi_save((post_j-1)/nthin) = c_xi_samp;
        }

        if (learn_c_tau & (mod_type == "triple")) {
          c_tau_save((post_j-1)/nthin) = c_tau_samp;
        }


        if (sv == false) {
          C0_save((post_j-1)/nthin) = C0_samp;
        } else {
          sv_mu_save((post_j-1)/nthin) = sv_para(0);
          sv_phi_save((post_j-1)/nthin) = sv_para(1);
          sv_sigma2_save((post_j-1)/nthin) = sv_para(2);
        }
      }

      // Conditionally store MH statistics
      if (learn_a_xi & bool(adaptive(0)) & (batch_pos(0) == (batch_sizes(0) - 2))){
        a_xi_sd_save(batch_nrs(0) - 1) = curr_sds(0);
        a_xi_acc_rate_save(batch_nrs(0) - 1) = arma::accu(batches.col(0))/batch_sizes(0);
      }
      if (learn_a_tau & bool(adaptive(1)) & (batch_pos(1) == (batch_sizes(1) - 2))){
        a_tau_sd_save(batch_nrs(1) - 1) = curr_sds(1);
        a_tau_acc_rate_save(batch_nrs(1) - 1) = arma::accu(batches.col(1))/batch_sizes(1);
      }
      if (learn_c_xi & bool(adaptive(2)) & (batch_pos(2) == (batch_sizes(2) - 2))){
        c_xi_sd_save(batch_nrs(2) - 1) = curr_sds(2);
        c_xi_acc_rate_save(batch_nrs(2) - 1) = arma::accu(batches.col(2))/batch_sizes(2);
      }
      if (learn_c_tau & bool(adaptive(3)) & (batch_pos(3) == (batch_sizes(3) - 2))){
        c_tau_sd_save(batch_nrs(3) - 1) = curr_sds(3);
        c_tau_acc_rate_save(batch_nrs(3) - 1) = arma::accu(batches.col(3))/batch_sizes(3);
      }
    }

    // Random sign switch
    if (niter > 1) {
      for (int i = 0; i < d; i++) {
        if (R::runif(0,1) > 0.5) {
          theta_sr_samp(i) = -theta_sr_samp(i);
        }
      }
    }

    // Increment progress bar
    if (arma::any(prog_rep_points == j)) {
      p.increment();
    }

    // Check for user interrupts
    if (j % 500 == 0) {
      Rcpp::checkUserInterrupt();
    }

    // Break loop if succesful is false
    if (!succesful) {
      break;
    }
  }

  // if niter == 1 && nburn == 0 we know that shrinkTVP is being used in the one step update mode - hence the return value changes
  // return everything as a nested list (due to size restrictions on Rcpp::Lists)
  if (niter == 1 && nburn == 0) {

    if(mod_type == "triple") {
      tau2_samp = tau2_til_samp;
      xi2_samp = xi2_til_samp;
    }

    return List::create(_["beta"] = (beta_nc_samp.each_col() % theta_sr_samp).each_col() + beta_mean_samp,
                        _["beta_mean_st"] = beta_mean_samp,
                        _["theta_sr_st"] = theta_sr_samp,
                        _["tau2_st"] = tau2_samp,
                        _["xi2_st"] = xi2_samp,
                        _["lambda2_st"] = lambda2_til_samp,
                        _["kappa2_st"] = kappa2_til_samp,
                        _["a_xi_st"] = a_xi_samp,
                        _["c_xi_st"] = c_xi_samp,
                        _["a_tau_st"] = a_tau_samp,
                        _["c_tau_st"] = c_tau_samp,
                        _["lambda2_B_st"] = lambda2_B_samp,
                        _["kappa2_B_st"] = kappa2_B_samp,
                        _["sigma2_st"] = sigma2_samp,
                        _["C0_st"] = C0_samp,
                        _["sv_mu_st"] = sv_para(0),
                        _["sv_phi_st"] = sv_para(1),
                        _["sv_sigma2_st"] = sv_para(2),
                        _["h0_st"] = h0_samp,
                        _["internals"] = List::create(
                          _["m_N"] = m_N_samp,
                          _["chol_C_N_inv"] = chol_C_N_inv_samp,
                          _["batch_pos_st"] = batch_pos,
                          _["batches_st"] = batches,
                          _["curr_sds_st"] = curr_sds,
                          _["batch_nrs_st"] = batch_nrs,
                          _["success"] = succesful,
                          _["fail"] = fail));
  } else {
    return List::create(_["beta"] = beta_save,
                        _["beta_mean"] = beta_mean_save.t(),
                        _["theta_sr"] = theta_sr_save.t(),
                        _["tau2"] = tau2_save.t(),
                        _["xi2"] = xi2_save.t(),
                        _["lambda2"] = lambda2_save.t(),
                        _["kappa2"] = kappa2_save.t(),
                        _["a_xi"] = a_xi_save,
                        _["c_xi"] = c_xi_save,
                        _["a_tau"] = a_tau_save,
                        _["c_tau"] = c_tau_save,
                        _["kappa2_B"] = kappa2_B_save,
                        _["lambda2_B"] = lambda2_B_save,
                        _["sigma2"] = sigma2_save,
                        _["C0"] = C0_save,
                        _["sv_mu"] = sv_mu_save,
                        _["sv_phi"] = sv_phi_save,
                        _["sv_sigma2"] = sv_sigma2_save,
                        _["MH_diag"] = List::create(
                          _["a_xi_sds"] = a_xi_sd_save,
                          _["a_xi_acc_rate"] = a_xi_acc_rate_save,
                          _["a_tau_sds"] = a_tau_sd_save,
                          _["a_tau_acc_rate"] = a_tau_acc_rate_save,
                          _["c_xi_sds"] = c_xi_sd_save,
                          _["c_xi_acc_rate"] = c_xi_acc_rate_save,
                          _["c_tau_sds"] = c_tau_sd_save,
                          _["c_tau_acc_rate"] = c_tau_acc_rate_save
                        ),
                        _["internals"] = List::create(
                          _["m_N"] = m_N_save,
                          _["chol_C_N_inv"] = chol_C_N_inv_save,
                          _["success_vals"] = List::create(
                            _["success"] = succesful,
                            _["fail"] = fail,
                            _["fail_iter"] = fail_iter)));
  }
}



