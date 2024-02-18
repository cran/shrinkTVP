#' Markov Chain Monte Carlo (MCMC) for time-varying parameter models with dynamic shrinkage
#'
#' \code{shrinkTVP} samples from the joint posterior distribution of the parameters of a time-varying
#' parameter model with dynamic triple gamma shrinkage, potentially including stochastic volatility (SV),
#' and returns the MCMC draws.
#'
#' For details concerning the algorithms please refer to the papers by Bitto and Frühwirth-Schnatter (2019),
#' Cadonna et al. (2020) and Knaus and Frühwirth-Schnatter (2023).
#' For more details on the package and the usage of the functions, see Knaus et al. (2021).
#'
#' @param formula object of class "formula": a symbolic representation of the model, as in the
#' function \code{lm}. For details, see \code{\link{formula}}.
#' @param data \emph{optional} data frame containing the response variable and the covariates. If not found in \code{data},
#' the variables are taken from \code{environment(formula)}, typically the environment from which \code{shrinkTVP}
#' is called. No \code{NA}s are allowed in the response variable and the covariates.
#' @param mod_type character string that reads either \code{"triple"}, \code{"double"} or \code{"ridge"}.
#' Determines whether the triple gamma, double gamma or ridge prior are used for \code{theta_sr} and \code{beta_mean}.
#' The default is "double".
#' @param niter positive integer, indicating the number of MCMC iterations
#' to perform, including the burn-in. Has to be larger than or equal to \code{nburn} + 2. The default value is 10000.
#' @param nburn non-negative integer, indicating the number of iterations discarded
#' as burn-in. Has to be smaller than or equal to \code{niter} - 2. The default value is \code{round(niter / 2)}.
#' @param nthin positive integer, indicating the degree of thinning to be performed. Every \code{nthin} draw is kept and returned.
#' The default value is 1, implying that every draw is kept.
#' @param learn_a_xi logical value indicating whether to learn a_xi, the spike parameter of the state variances.
#' Ignored if \code{mod_type} is set to \code{"ridge"}. The default value is \code{TRUE}.
#' @param learn_a_tau logical value indicating whether to learn a_tau, the spike parameter of the mean of the
#' initial values of the states. Ignored if \code{mod_type} is set to \code{"ridge"}. The default value is \code{TRUE}.
#' @param a_xi positive, real number, indicating the (fixed) value for a_xi. Ignored if
#' \code{learn_a_xi} is \code{TRUE} or \code{mod_type} is set to \code{"ridge"}. The default value is 0.1.
#' @param a_tau positive, real number, indicating the (fixed) value for a_tau. Ignored if
#' \code{learn_a_tau} is \code{TRUE} or \code{mod_type} is set to \code{"ridge"}. The default value is 0.1.
#' @param learn_c_xi logical value indicating whether to learn c_xi, the tail parameter of the state variances.
#' Ignored if \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_xi} is set to \code{TRUE}.
#' The default value is \code{TRUE}.
#' @param learn_c_tau logical value indicating whether to learn c_tau, the tail parameter of the mean of the
#' initial values of the states. Ignored if \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_tau} is set to \code{TRUE}.
#' The default value is \code{TRUE}.
#' @param c_xi positive, real number, indicating the (fixed) value for c_xi. Ignored if
#' \code{learn_c_xi} is \code{TRUE}, \code{mod_type} is not set to \code{"triple"} or \code{a_eq_c_xi} is set to \code{TRUE}.
#' The default value is 0.1.
#' @param c_tau positive, real number, indicating the (fixed) value for c_tau. Ignored if
#' \code{learn_c_xi} is \code{TRUE}, \code{mod_type} is not set to \code{"triple"}  or \code{a_eq_c_tau} is set to \code{TRUE}.
#' The default value is 0.1.
#' @param a_eq_c_xi logical value indicating whether to force \code{a_xi} and \code{c_xi} to be equal.
#' If set to \code{TRUE}, \code{beta_a_xi} and \code{alpha_a_xi} are used as the hyperparameters and \code{beta_c_xi} and \code{alpha_c_xi} are ignored.
#' Ignored if \code{mod_type} is not set to \code{"triple"}. The default value is \code{FALSE}.
#' @param a_eq_c_tau logical value indicating whether to force \code{a_tau} and \code{c_tau} to be equal.
#' If set to \code{TRUE}, \code{beta_a_tau} and \code{alpha_a_tau} are used as the hyperparameters and \code{beta_c_tau} and \code{alpha_c_tau} are ignored.
#' Ignored if \code{mod_type} is not set to \code{"triple"}. The default value is \code{FALSE}.
#' @param learn_kappa2_B logical value indicating whether to learn kappa2_B, the global level of shrinkage for
#' the state variances. The default value is \code{TRUE}.
#' @param learn_lambda2_B logical value indicating whether to learn the lambda2_B parameter,
#' the global level of shrinkage for the mean of the initial values of the states. The default value is \code{TRUE}.
#' @param kappa2_B positive, real number, indicating the (fixed) value for kappa2_B. Ignored if
#' \code{learn_kappa2_B} is \code{TRUE}. The default value is 20.
#' @param lambda2_B positive, real number, indicating the (fixed) value for lambda2_B. Ignored if
#' \code{learn_lambda2_B} is \code{TRUE}. The default value is 20.
#' @param a_psi positive, real number, or a vector of length equal to the number of covariates containing
#' positive, real numbers. Indicates the value for a_psi, which is the pole parameter of the dynamic triple gamma.
#' The default value is 0.5.
#' @param c_psi positive, real number, or a vector of length equal to the number of covariates containing
#' positive, real numbers. Indicates the value for c_psi, which is the tail parameter of the dynamic triple gamma.
#' The default value is 0.5.
#' @param iid logical value indicating whether the innovations are assumed to be independent and identically distributed.
#' If set to \code{TRUE}, the innovations are assumed to be a priori iid triple gamma. If set to \code{FALSE},
#' the prior on the innovations is the dynamic triple gamma specification of Knaus and Frühwirth-Schnatter (2023).
#' The default value is \code{FALSE}.
#' @param shrink_inter logical value indicating whether to dynamically shrink the intercept. Note that shrinkage is still applied
#' to the \code{theta_sr} and \code{beta_mean} associated with the intercept. The intercept column is automatically determined
#' by the function and does not have to be included in the formula. This is done by finding the column that contains only 1s.
#' The default value is \code{TRUE}.
#' @param hyperprior_param \emph{optional} named list containing hyperparameter values.
#' Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown.
#' All hyperparameter values have to be positive, real numbers. The following hyperparameters can be
#' supplied:
#' \itemize{
#' \item \code{c0}: The default value is 2.5.
#' \item \code{g0}: The default value is 5.
#' \item \code{G0}: The default value is 5 / (2.5 - 1).
#' \item \code{e1}: The default value is 0.001.
#' \item \code{e2}: The default value is 0.001.
#' \item \code{d1}: The default value is 0.001.
#' \item \code{d2}: The default value is 0.001.
#' \item \code{alpha_a_xi}: The default value is 5.
#' \item \code{alpha_a_tau}: The default value is 5.
#' \item \code{beta_a_xi}: The default value is 10.
#' \item \code{beta_a_tau}: The default value is 10.
#' \item \code{alpha_c_xi}: The default value is 5.
#' \item \code{alpha_c_tau}: The default value is 5.
#' \item \code{beta_c_xi}: The default value is 2.
#' \item \code{beta_c_tau}: The default value is 2.
#' \item \code{a_rho}: The default value is 2.
#' \item \code{b_rho}: The default value is 0.95.
#' \item \code{alpha_rho}: The default value is 0.5.
#' \item \code{beta_rho}: The default value is 3.
#' }
#' @param display_progress logical value indicating whether the progress bar and other informative output should be
#' displayed. The default value is \code{TRUE}.
#' @param sv logical value indicating whether to use stochastic volatility for the error of the observation
#' equation. For details please see \code{\link{stochvol}}, in particular \code{\link{svsample}}. The default value is
#' \code{FALSE}.
#' @param sv_param \emph{optional} named list containing hyperparameter values for the stochastic volatility
#' parameters. Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown. Ignored if
#' \code{sv} is \code{FALSE}. The following elements can be supplied:
#' \itemize{
#' \item \code{Bsigma_sv}: positive, real number. The default value is 1.
#' \item \code{a0_sv}: positive, real number. The default value is 5.
#' \item \code{b0_sv}: positive, real number. The default value is 1.5.
#' \item \code{bmu}: real number. The default value is 0.
#' \item \code{Bmu}: real number. larger than 0. The default value is 1.
#' }
#' @param MH_tuning \emph{optional} named list containing values used to tune the MH steps for \code{a_xi}, \code{a_tau},
#' \code{c_xi}, and \code{c_tau}. Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown.
#' The arguments for \code{a_xi}(\code{a_tau}) are only used if \code{learn_a_xi}(\code{learn_a_tau})
#' is set to \code{TRUE} and \code{mod_type} is not equal to \code{"ridge"}. The arguments for \code{c_xi}(\code{c_tau}) are only
#' used if \code{learn_c_xi}(\code{learn_c_tau}) is set to \code{TRUE} and \code{mod_type} is equal to \code{"triple"}. Arguments ending in "adaptive" are
#' logical values indicating whether or not to make the MH step for the respective parameter adaptive. Arguments ending in "tuning_par" serve two different purposes.
#' If the respective MH step is not set to be adaptive, it acts as the standard deviation of the proposal distribution. If the respective MH step
#' is set to be adaptive, it acts as the initial standard deviation. Arguments ending in "target_rate" define the acceptance rate the algorithm aims to achieve.
#' Arguments ending in "max_adapt" set the maximum value by which the logarithm of the standard deviation of the proposal distribution is adjusted. Finally,
#' arguments ending in "batch_size" set the batch size after which the standard deviation of the proposal distribution is adjusted.
#' The following elements can be supplied:
#' \itemize{
#' \item \code{a_xi_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{a_xi_tuning_par}: positive, real number. The default value is 1.
#' \item \code{a_xi_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{a_xi_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{a_xi_batch_size}: positive integer. The default value is 50.
#' \item \code{a_tau_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{a_tau_tuning_par}: positive, real number. The default value is 1.
#' \item \code{a_tau_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{a_tau_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{a_tau_batch_size}: positive integer. The default value is 50.
#' \item \code{c_xi_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{c_xi_tuning_par}: positive, real number. The default value is 1.
#' \item \code{c_xi_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{c_xi_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{c_xi_batch_size}: positive integer. The default value is 50.
#' \item \code{c_tau_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{c_tau_tuning_par}: positive, real number. The default value is 1.
#' \item \code{c_tau_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{c_tau_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{c_tau_batch_size}: positive integer. The default value is 50.
#' \item \code{rho_adaptive}: logical value. The default is \code{TRUE}.
#' \item \code{rho_tuning_par}: positive, real number. The default value is 1.
#' \item \code{rho_target_rate}: positive, real number, between 0 and 1. The default value is 0.44.
#' \item \code{rho_max_adapt}: positive, real number. The default value is 0.01.
#' \item \code{rho_batch_size}: positive integer. The default value is 50.
#' }
#' @param starting_vals \emph{optional} named list containing the values at which the MCMC algorithm will be initialized. In the
#' following \code{d} refers to the number of covariates, including the intercept and expanded factors.
#' Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown. The following elements can be supplied:
#' \itemize{
#' \item \code{beta_mean_st}: vector of length \code{d} containing single numbers. The default is \code{rep(0, d)}.
#' \item \code{theta_sr_st}: vector of length \code{d} containing single, positive numbers. The default is \code{rep(1, d)}.
#' \item \code{tau2_st}: vector of length \code{d} containing single, positive numbers. The default is \code{rep(1, d)}.
#' \item \code{xi2_st}: vector of length \code{d} containing single, positive numbers. The default is \code{rep(1, d)}.
#' \item \code{kappa2_st}: vector of length \code{d} containing single, positive numbers. The default is \code{rep(1, d)}.
#' \item \code{lambda2_st}: vector of length \code{d} containing single, positive numbers. The default is \code{rep(1, d)}.
#' \item \code{kappa2_B_st}: positive, real number. The default value is 20.
#' \item \code{lambda2_B_st}: positive, real number. The default value is 20.
#' \item \code{a_xi_st}: positive, real number. The default value is 0.1.
#' \item \code{a_tau_st}: positive, real number. The default value is 0.1.
#' \item \code{c_xi_st}: positive, real number. The default value is 0.1. Note that the prior for \code{c_xi} is restricted to (0, 0.5).
#' \item \code{c_tau_st}: positive, real number. The default value is 0.1. Note that the prior for \code{c_tau} is restricted to (0, 0.5).
#' \item \code{sv_mu_st}: real number. The default value is -10.
#' \item \code{sv_phi_st}: positive, real number between -1 and 1. The default value is 0.5.
#' \item \code{sv_sigma2_st }: positive, real number. The default value is 1.
#' \item \code{C0_st}: positive, real number. The default value is 1.
#' \item \code{sigma2_st}: positive, real number if \code{sv} is \code{FALSE}, otherwise a vector of positive, real numbers of length \code{N}. The default value is 1 or a vector thereof.
#' \item \code{h0_st}: real number. The default value is 0.
#' \item \code{lambda_0_st} vector of length \code{d} containing positive, real numbers. The default value is \code{rep(1, d)}.
#' \item \code{rho_st}: vector of length \code{d} containing real numbers between 0 and \code{b_rho}. The default value is \code{rep(max(0.1, hyperprior_param$b_rho - 0.1), d)}.
#' }
#' @return The value returned is a list object of class \code{shrinkTVP} containing
#' \item{\code{beta}}{\code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior distribution of the centered
#' states, one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.tvp} object.}
#' \item{\code{beta_mean}}{\code{mcmc} object containing the parameter draws from the posterior distribution of beta_mean.}
#' \item{\code{theta_sr}}{\code{mcmc} object containing the parameter draws from the posterior distribution of the square root of theta.}
#' \item{\code{tau2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of tau2.}
#' \item{\code{xi2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of xi2.}
#' \item{\code{psi}}{\code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior distribution of psi,
#' one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.tvp} object.}
#' \item{\code{lambda_p}}{\code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior distribution of lambda_p,
#' one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.tvp} object.}
#' \item{\code{kappa_p}}{\emph{(optional)} \code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior
#'  distribution of kappa_p, one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.tvp} object.
#' Not returned if \code{iid} is not \code{TRUE}.}
#' \item{\code{lambda2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of lambda2.
#' Not returned if \code{mod_type} is not \code{"triple"}.}
#' \item{\code{kappa2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of kappa2.
#' Not returned if \code{mod_type} is not \code{"triple"}.}
#' \item{\code{a_xi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_xi.
#' Not returned if \code{learn_a_xi} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{a_tau}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_tau.
#' Not returned if \code{learn_a_tau} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{c_xi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of c_xi.
#' Not returned if \code{learn_c_xi} is \code{FALSE} or \code{mod_type} is not \code{"triple"}.}
#' \item{\code{c_tau}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of c_tau.
#' Not returned if \code{learn_c_tau} is \code{FALSE} or \code{mod_type} is not \code{"triple"}.}
#' \item{\code{lambda2_B}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of lambda2_B.
#' Not returned if \code{learn_lambda2_B} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{kappa2_B}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of kappa2_B.
#' Not returned if \code{learn_kappa2_B} is \code{FALSE} or \code{mod_type} is \code{"ridge"}.}
#' \item{\code{rho}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of rho.
#' Not returned if \code{iid} is not \code{TRUE}.}
#' \item{\code{sigma2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of \code{sigma2}.
#' If \code{sv} is \code{TRUE}, \code{sigma2} is additionally an \code{mcmc.tvp} object.}
#' \item{\code{C0}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of C0.
#' Not returned if \code{sv} is \code{TRUE}.}
#' \item{\code{sv_mu}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the mu
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{sv_phi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the phi
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{sv_sigma2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the sigma2
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{MH_diag}}{\emph{(optional)} named list containing statistics for assessing MH performance. Not returned if no MH steps are required
#' or none of them are specified to be adaptive.}
#' \item{\code{internals}}{\code{list} object containing two arrays that are required for calculating the LPDS.}
#' \item{\code{priorvals}}{\code{list} object containing hyperparameter values of the prior distributions, as specified by the user.}
#' \item{\code{model}}{\code{list} object containing the model matrix, model response and formula used.}
#' \item{\code{summaries}}{\code{list} object containing a collection of summary statistics of the posterior draws.}
#'
#' To display the output, use \code{plot} and \code{summary}. The \code{summary} method displays the specified prior values stored in
#' \code{priorvals} and the posterior summaries stored in \code{summaries}, while the \code{plot} method calls \code{coda}'s \code{plot.mcmc}
#' or the \code{plot.mcmc.tvp} method. Furthermore, all functions that can be applied to \code{coda::mcmc} objects
#' (e.g. \code{coda::acfplot}) can be applied to all output elements that are \code{coda} compatible.
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#'
#' @seealso \code{\link{plot.shrinkTVP}}, \code{\link{plot.mcmc.tvp}}
#' @references Bitto, A., & Frühwirth-Schnatter, S. (2019). "Achieving shrinkage in a time-varying parameter model framework."
#' \emph{Journal of Econometrics}, 210(1), 75-97. <doi:10.1016/j.jeconom.2018.11.006>
#'
#' Cadonna, A., Frühwirth-Schnatter, S., & Knaus, P. (2020). "Triple the Gamma—A Unifying Shrinkage Prior for Variance and Variable Selection in Sparse State Space and TVP Models."
#' \emph{Econometrics}, 8(2), 20. <doi:10.3390/econometrics8020020>
#'
#' Knaus, P., Bitto-Nemling, A., Cadonna, A., & Frühwirth-Schnatter, S. (2021) "Shrinkage in the Time-Varying Parameter Model Framework Using the \code{R} Package \code{shrinkTVP}."
#' \emph{Journal of Statistical Software} 100(13), 1–32. <doi:10.18637/jss.v100.i13>
#'
#' Knaus, P., & Frühwirth-Schnatter, S. (2023). "The Dynamic Triple Gamma Prior as a Shrinkage Process Prior for Time-Varying Parameter Models." arXiv preprint arXiv:2312.10487. <doi:10.48550/arXiv.2312.10487>
#' @examples
#' \donttest{
#'
## Example 1, stock model
#'set.seed(123)
#'sim <- simTVP(DTG = TRUE, theta = c(0, 1, 0), beta_mean = c(1, 1, 0), rho = 0.95, c_psi = 2)
#'data <- sim$data
#'
#'## Example 1, match the true underlying process
#'res <- shrinkDTVP(y ~ x1 + x2, data = data, c_psi = 2)
#'# summarize output
#'summary(res)
#'
#'## Example 2, dynamic horseshoe
#'res <- shrinkDTVP(y ~ x1 + x2, data = data)
#'
#'
#'## Example 3, modify hyperparameters
#'res <- shrinkDTVP(y ~ x1 + x2, data = data,
#'                  hyperprior_param = list(a_rho = 1,
#'                                          alpha_rho = 0.5,
#'                                          beta_rho = 0.5))
#'
#' }
#'
#' @export
shrinkDTVP <- function(formula,
                       data,
                       mod_type = "double",
                       niter = 10000,
                       nburn = round(niter / 2),
                       nthin = 1,
                       learn_a_xi = TRUE,
                       learn_a_tau = TRUE,
                       a_xi = 0.1,
                       a_tau = 0.1,
                       learn_c_xi = TRUE,
                       learn_c_tau = TRUE,
                       c_xi = 0.1,
                       c_tau = 0.1,
                       a_eq_c_xi = FALSE,
                       a_eq_c_tau = FALSE,
                       learn_kappa2_B = TRUE,
                       learn_lambda2_B = TRUE,
                       kappa2_B = 20,
                       lambda2_B = 20,
                       a_psi = 0.5,
                       c_psi = 0.5,
                       iid = FALSE,
                       shrink_inter = TRUE,
                       hyperprior_param,
                       display_progress = TRUE,
                       sv = FALSE,
                       sv_param,
                       MH_tuning,
                       starting_vals){


  # Input checking ----------------------------------------------------------

  # check that mod_type was specified correctly
  if (!mod_type %in% c("triple", "double", "ridge")) {
    stop("mod_type has to be a string equal to 'triple', 'double' or 'ridge'")
  }

  # default hyperparameter values
  default_hyper <- list(c0 = 2.5,
                        g0 = 5,
                        G0 = 5 / (2.5 - 1),
                        e1 = 0.001,
                        e2 = 0.001,
                        d1 = 0.001,
                        d2 = 0.001,
                        beta_a_xi = 10,
                        beta_a_tau = 10,
                        alpha_a_xi = 5,
                        alpha_a_tau = 5,
                        beta_c_xi = 2,
                        beta_c_tau = 2,
                        alpha_c_xi = 5,
                        alpha_c_tau = 5,
                        a_rho = 2,
                        b_rho = .95,
                        alpha_rho = .5,
                        beta_rho = 3)

  # default sv params
  default_hyper_sv <- list(Bsigma_sv = 1,
                           a0_sv = 5,
                           b0_sv = 1.5,
                           bmu = 0,
                           Bmu = 1)

  # default tuning parameters
  default_tuning_par <- list(a_xi_adaptive = TRUE,
                             a_xi_tuning_par = 1,
                             a_xi_target_rate = 0.44,
                             a_xi_max_adapt = 0.01,
                             a_xi_batch_size = 50,
                             a_tau_adaptive = TRUE,
                             a_tau_tuning_par = 1,
                             a_tau_target_rate = 0.44,
                             a_tau_max_adapt = 0.01,
                             a_tau_batch_size = 50,
                             c_xi_adaptive = TRUE,
                             c_xi_tuning_par = 1,
                             c_xi_target_rate = 0.44,
                             c_xi_max_adapt = 0.01,
                             c_xi_batch_size = 50,
                             c_tau_adaptive = TRUE,
                             c_tau_tuning_par = 1,
                             c_tau_target_rate = 0.44,
                             c_tau_max_adapt = 0.01,
                             c_tau_batch_size = 50,
                             rho_adaptive = TRUE,
                             rho_tuning_par = 1,
                             rho_target_rate = 0.44,
                             rho_max_adapt = 0.01,
                             rho_batch_size = 50)

  # Change tuning parameter values if user overwrites them
  if (missing(MH_tuning)){
    MH_tuning <- default_tuning_par
  } else {
    MH_tuning <- list_merger(default_tuning_par, MH_tuning)
  }


  # Change hyperprior values if user overwrites them
  if (missing(hyperprior_param)){
    hyperprior_param <- default_hyper
  } else {
    hyperprior_param <- list_merger(default_hyper, hyperprior_param)
  }


  # Same procedure for sv_param
  if (missing(sv_param) | sv == FALSE){
    sv_param <- default_hyper_sv
  } else {
    sv_param <- list_merger(default_hyper_sv, sv_param)
  }

  # Check if all numeric inputs are correct
  to_test_num <- list(lambda2_B = lambda2_B,
                      kappa2_B = kappa2_B,
                      a_xi = a_xi,
                      a_tau = a_tau,
                      c_xi = c_xi,
                      c_tau = c_tau)


  if (missing(hyperprior_param) == FALSE){
    to_test_num <- c(to_test_num, hyperprior_param)
  }

  if (missing(sv_param) == FALSE){
    to_test_num <- c(to_test_num, sv_param[names(sv_param) != "bmu"])
  }

  if (missing(MH_tuning) == FALSE){
    to_test_num <- c(to_test_num, MH_tuning[!grepl("(batch|adaptive)", names(MH_tuning))])
  }

  bad_inp <- sapply(to_test_num, numeric_input_bad)


  if (any(bad_inp)){
    stand_names <- c(names(default_hyper), names(default_hyper_sv), "lambda2_B", "kappa2_B", "a_xi", "a_tau", "c_xi", "c_tau")
    bad_inp_names <- names(to_test_num)[bad_inp]
    bad_inp_names <- bad_inp_names[bad_inp_names %in% stand_names]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a real, positive number"))
  }

  # Check bmu seperately
  if (!is.numeric(sv_param$bmu) | !is.scalar(sv_param$bmu)){
    stop("bmu has to be a real number")
  }

  # Check the adapt_rates seperately
  if (any(0 > MH_tuning[grepl("rate", names(MH_tuning))] | MH_tuning[grepl("rate", names(MH_tuning))] > 1)) {
    stop("all target_rate parameters in MH_tuning have to be > 0 and < 1")
  }

  # Check if all integer inputs are correct
  to_test_int <- c(niter = niter,
                   nburn = nburn,
                   nthin = nthin,
                   #p = p,
                   MH_tuning[grepl("batch", names(MH_tuning))])

  bad_int_inp <- sapply(to_test_int, int_input_bad)

  if (any(bad_int_inp)){
    bad_inp_names <- names(to_test_int)[bad_int_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive integer"))

  }

  if ((niter - nburn) < 2){
    stop("niter has to be larger than or equal to nburn + 2")
  }

  if (nthin == 0){
    stop("nthin can not be 0")
  }

  if ((niter - nburn)/2 < nthin){
    stop("nthin can not be larger than (niter - nburn)/2")
  }

  # Check if all boolean inputs are correct
  to_test_bool <- c(learn_lambda2_B = learn_lambda2_B,
                    learn_kappa2_B = learn_kappa2_B,
                    learn_a_xi = learn_a_xi,
                    learn_a_tau = learn_a_tau,
                    display_progress = display_progress,
                    sv = sv,
                    MH_tuning[grepl("adaptive", names(MH_tuning))],
                    iid,
                    shrink_inter)

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- names(to_test_bool)[bad_bool_inp]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single logical value"))

  }

  # Check if formula is a formula
  if (inherits(formula, "formula") == FALSE){
    stop("formula is not of class formula")
  }




  # Formula interface -------------------------------------------------------


  mf <- match.call(expand.dots = FALSE)
  m <- match(x = c("formula", "data"), table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  # Create Vector y
  y <- model.response(mf, "numeric")
  mt <- attr(x = mf, which = "terms")
  # Create Matrix X with dummies and transformations
  x <- model.matrix(object = mt, data = mf)

  # Check that there are no NAs in y and x
  if (any(is.na(y))) {
    stop("No NA values are allowed in response variable")
  }

  if (any(is.na(x))){
    stop("No NA values are allowed in covariates")
  }

  # Get the index
  if (missing(data)){
    index <- zoo::index(y)
  } else {
    index <- zoo::index(data)
  }

  p <- 0
  # Create lagged values of y and add to list of regressors
  if (p != 0){

    # Create new x matrix with added in lagged values
    x <- cbind(x, mlag(y, p))[(p + 1):nrow(x), ]

    # Rename newly added columns to ar1...arp
    colnames(x)[(ncol(x) - p + 1):ncol(x)] <- paste0("ar", 1:p)

    # Cut y length to appropriate length
    y <- y[(p + 1):length(y)]

    index <- index[(p + 1):length(index)]

  }

  colnames(x)[colnames(x) == "(Intercept)"] <- "Intercept"
  d <- dim(x)[2]

  # Find intercept column, if it exists
  # Note that we do not simply check the first column, or the column named "Intercept"
  # as the user may have added an intercept column manually

  inter_column <- which(apply(x, 2, function(x) all(x == 1)))

  if (length(inter_column) > 1) {
    stop("More than one intercept column (i.e. column with all elements equal to 1) found")
  }

  if (shrink_inter) {
    if (length(inter_column) == 0) {
      warning("No intercept column found, but shrink_inter is TRUE", immediate. = TRUE)
    }
    inter_column <- NA_integer_
  } else {
    if (length(inter_column) == 0) {
      inter_column <- NA_integer_
    }
  }


  # Fuse user starting vals with standard ones
  default_starting_vals <- list(psi_st = matrix(1, d, length(y)),
                                beta_mean_st = rep(0, d),
                                theta_sr_st = rep(1, d),
                                tau2_st = rep(1, d),
                                xi2_st = rep(1, d),
                                kappa2_st = rep(1, d),
                                lambda2_st = rep(1, d),
                                kappa2_B_st = 20,
                                lambda2_B_st = 20,
                                a_xi_st = 0.1,
                                a_tau_st = 0.1,
                                c_xi_st = 0.1,
                                c_tau_st = 0.1,
                                sv_mu_st = -10,
                                sv_phi_st = 0.5,
                                sv_sigma2_st = 1,
                                C0_st = 1,
                                sigma2_st = 1,
                                h0_st = 0,
                                rho_st = rep(max(0.1, hyperprior_param$b_rho - 0.1), d),
                                lambda_0_st = rep(1, d))

  if (sv == TRUE){
    default_starting_vals$sigma2_st <- rep(1, length(y))
  }

  hyperprior_param$alpha_rho <- rep(hyperprior_param$alpha_rho, d)
  hyperprior_param$beta_rho <- rep(hyperprior_param$beta_rho, d)

  # Change starting values of MCMC algorithm if user overwrites them
  if (missing(starting_vals)){
    starting_vals <- default_starting_vals
  } else {
    starting_vals <- list_merger(default_starting_vals, starting_vals)
  }

  # Input check starting vals

  # Check matrix valued inputs
  # This is a bit over-engineered, to allow for more matrix valued inputs in the future
  matrix_valued <- c("psi_st")

  # First, check that they are matrices
  bad_matrix <- sapply(starting_vals[matrix_valued], function(x) !is.matrix(x))
  if (any(bad_matrix)){
    bad_matrix_names <- matrix_valued[bad_matrix]
    stop(paste0(paste(bad_matrix_names, collapse = ", "),
                ifelse(length(bad_matrix_names) == 1, " has", " have"),
                " to be of class matrix"))
  }

  # Next, check that they have the right dimensions
  bad_matrix_dim <- sapply(starting_vals[matrix_valued], function(x) any(dim(x) != c(d, length(y))))
  if (any(bad_matrix_dim)){
    bad_matrix_dim_names <- matrix_valued[bad_matrix_dim]
    stop(paste0(paste(bad_matrix_dim_names, collapse = ", "),
                ifelse(length(bad_matrix_dim_names) == 1, " has", " have"),
                " to be of dimension d = ", d, " x N = ", length(y)))
  }

  # Check length of vectors
  vec_valued <- c("beta_mean_st",
                  "theta_sr_st",
                  "tau2_st",
                  "xi2_st",
                  "kappa2_st",
                  "lambda2_st",
                  "lambda_0_st",
                  "rho_st")

  bad_length <- sapply(starting_vals[vec_valued], function(x) length(x) != d)

  if (any(bad_length)){
    bad_length_names <- vec_valued[bad_length]
    stop(paste0(paste(bad_length_names, collapse = ", "),
                ifelse(length(bad_length_names) == 1, " has", " have"),
                " to be of length ", d))

  }

  # check sigma2_st seperately
  if (sv == TRUE) {
    if (length(starting_vals$sigma2_st) != length(y)) {
      stop("sigma2_st has to be the same length as y if sv is TRUE")
    }

    num_input_bad <- sapply(starting_vals$sigma2_st, numeric_input_bad)

    if (any(num_input_bad)) {
      stop("sigma2_st may only contain real, positive numbers")
    }
  } else if (numeric_input_bad(starting_vals$sigma2_st)) {
    stop("sigma2_st has to be a real, positive number")
  }

  # Check content of vectors
  vec_valued_pos <- vec_valued[vec_valued != "beta_mean_st"]
  bad_content <- sapply(starting_vals[vec_valued_pos], function(x) any(sapply(x, numeric_input_bad)))

  if (any(bad_content)) {
    bad_content_names <- vec_valued_pos[bad_content]
    stop(paste0(paste(bad_content_names, collapse = ", "), " may only contain real, positive numbers"))

  }

  # Check beta_mean_st seperately
  if (any(sapply(starting_vals$beta_mean_st, numeric_input_bad_))) {
    stop("beta_mean_st may only contain real numbers")
  }

  # Check rho_st seperately
  if (any(sapply(starting_vals$rho_st, function(x) x < 0 | x > 1))) {
    stop("rho_st may only contain real numbers between 0 and 1")
  }

  # Check that rho_st is in the support of the prior (i.e. smaller than b_rho)
  if (any(sapply(starting_vals$rho_st, function(x) x > hyperprior_param$b_rho))) {
    stop("At least one element of rho_st is larger than b_rho, the upper bound of the support of the prior for rho")
  }

  # Check single values
  to_check_num_st <- names(starting_vals)[!names(starting_vals) %in% c(vec_valued,  matrix_valued, "h0_st", "sv_mu_st", "sigma2_st", "do_dat_G")]
  bad_num_st <- sapply(starting_vals[to_check_num_st], numeric_input_bad)

  if (any(bad_num_st)) {
    bad_num_names <- to_check_num_st[bad_num_st]
    stop(paste0(paste(bad_num_names, collapse = ", "),
                ifelse(length(bad_num_names) == 1, " has", " have"),
                " to be a real, positive number"))

  }

  # Check h0_st and sv_mu_st seperately
  if (numeric_input_bad_(starting_vals$h0_st)) {
    stop("h0_st has to be a real number")
  }

  if (numeric_input_bad_(starting_vals$sv_mu_st)) {
    stop("sv_mu_st has to be a real number")
  }

  # Check that sv_phi_st falls between -1 and 1
  if (abs(starting_vals$sv_phi_st) >= 1) {
    stop("sv_phi_st has to be between -1 and 1")
  }

  # Check a_psi and c_psi
  # First check length
  if (!length(a_psi) %in% c(1, d)) {
    stop(paste0("a_psi has to be of length 1 or d = ", d))
  }
  if (!length(c_psi) %in% c(1, d)) {
    stop(paste0("c_psi has to be of length 1 or d = ", d))
  }

  # Then check content
  if (any(sapply(a_psi, numeric_input_bad))) {
    stop("a_psi may only contain positive real numbers")
  }
  if (any(sapply(c_psi, numeric_input_bad))) {
    stop("c_psi may only contain positive real numbers")
  }


  # Warn user if c_xi/c_tau outside the range of (0, 0.5)
  if (starting_vals$c_xi_st >= 0.5) {
    warning("c_xi_st is >= 0.5, which means the algorithm will not move away from the starting value", immediate. = TRUE)
  }
  if (starting_vals$c_tau_st >= 0.5) {
    warning("c_tau_st is >= 0.5, which means the algorithm will not move away from the starting value", immediate. = TRUE)
  }


  # Run sampler -------------------------------------------------------------

  if (length(a_psi) == 1) {
    a_psi <- rep(a_psi, d)
  }

  if (length(c_psi) == 1) {
    c_psi <- rep(c_psi, d)
  }

  runtime <- system.time({
    suppressWarnings({
      res <- shrinkDTVP_cpp(y,
                            x,
                            mod_type,
                            iid,
                            niter,
                            nburn,
                            nthin,
                            hyperprior_param$c0,
                            hyperprior_param$g0,
                            hyperprior_param$G0,
                            hyperprior_param$d1,
                            hyperprior_param$d2,
                            hyperprior_param$e1,
                            hyperprior_param$e2,
                            learn_lambda2_B,
                            learn_kappa2_B,
                            lambda2_B,
                            kappa2_B,
                            learn_a_xi,
                            learn_a_tau,
                            a_xi,
                            a_tau,
                            learn_c_xi,
                            learn_c_tau,
                            c_xi,
                            c_tau,
                            a_eq_c_xi,
                            a_eq_c_tau,
                            MH_tuning$a_xi_tuning_par,
                            MH_tuning$a_tau_tuning_par,
                            MH_tuning$c_xi_tuning_par,
                            MH_tuning$c_tau_tuning_par,
                            hyperprior_param$beta_a_xi,
                            hyperprior_param$beta_a_tau,
                            hyperprior_param$alpha_a_xi,
                            hyperprior_param$alpha_a_tau,
                            hyperprior_param$beta_c_xi,
                            hyperprior_param$beta_c_tau,
                            hyperprior_param$alpha_c_xi,
                            hyperprior_param$alpha_c_tau,
                            hyperprior_param$alpha_rho,
                            hyperprior_param$beta_rho,
                            a_psi,
                            c_psi,
                            hyperprior_param$a_rho,
                            hyperprior_param$b_rho,
                            inter_column,
                            display_progress,
                            sv,
                            sv_param$Bsigma_sv,
                            sv_param$a0_sv,
                            sv_param$b0_sv,
                            sv_param$bmu,
                            sv_param$Bmu,
                            MH_tuning$rho_adaptive,
                            MH_tuning$rho_tuning_par,
                            MH_tuning$rho_target_rate,
                            MH_tuning$rho_max_adapt,
                            MH_tuning$rho_batch_size,
                            unlist(MH_tuning[grep("adaptive", names(MH_tuning))]),
                            unlist(MH_tuning[grep("target", names(MH_tuning))]),
                            unlist(MH_tuning[grep("max", names(MH_tuning))]),
                            unlist(MH_tuning[grep("size", names(MH_tuning))]),
                            starting_vals)
    })
  })


  # Throw an error if the sampler failed
  if (res$internals$success_vals$success == FALSE){
    stop(paste0("The sampler failed at iteration ",
                res$internals$success_vals$fail_iter,
                " while trying to ",
                res$internals$success_vals$fail, ". ",
                "Try rerunning the model. ",
                "If the sampler fails again, try changing the prior to be more informative. ",
                "If the problem still persists, please contact the maintainer: ",
                maintainer("shrinkTVP")))
  } else {
    res$internals$success_vals <- NULL
  }

  # Post process sampler results --------------------------------------------

  # Unpack list of lists used to circumvent cap on number of elements in a list in Rcpp
  res <- c(res[1:(which(names(res) == "to_unpack") - 1)],
           res[["to_unpack"]],
           res[(which(names(res) == "to_unpack") + 1):length(names(res))])

  if (display_progress == TRUE){
    cat("Timing (elapsed): ", file = stderr())
    cat(runtime["elapsed"], file = stderr())
    cat(" seconds.\n", file = stderr())
    cat(round( (niter + nburn) / runtime[3]), "iterations per second.\n\n", file = stderr())
    cat("Converting results to coda objects and summarizing draws... ", file = stderr())
  }

  # Collapse sigma2 to single vector if sv=FALSE
  if (sv == FALSE){
    res$sigma2 <- matrix(res$sigma2[1, 1, ], ncol = 1)
  }

  # Remove empty storage elements
  res[sapply(res, function(x) 0 %in% dim(x))] <- NULL
  res$MH_diag[sapply(res$MH_diag, function(x) 0 %in% dim(x))] <- NULL

  if (a_eq_c_tau == TRUE) {
    res$c_tau <- NULL
  }

  if (a_eq_c_xi == TRUE) {
    res$c_xi <- NULL
  }


  # Create object to hold prior values
  res$priorvals <- c(hyperprior_param,
                     sv_param,
                     a_xi = a_xi,
                     a_tau = a_tau,
                     c_xi = c_xi,
                     c_tau = c_tau,
                     lambda2_B = lambda2_B,
                     kappa2_B = kappa2_B)

  # Add data to output
  res[["model"]] <- list()
  res$model$x <- x
  res$model$y <- y
  res$model$formula <- formula
  res$model$xlevels <- .getXlevels(mt, mf)
  res$model$terms <- mt

  res$summaries <- list()

  # add attributes to the individual objects if they are distributions or individual statistics
  nsave <- floor((niter - nburn)/nthin)
  for (i in names(res)){

    attr(res[[i]], "type") <- ifelse(nsave %in% dim(res[[i]]), "sample", "stat")

    # Name each individual sample for plotting frontend
    if (attr(res[[i]], "type") == "sample"){

      if (dim(res[[i]])[2] == d){

        colnames(res[[i]]) <- paste0(i, "_",  colnames(x))

      } else if (dim(res[[i]])[2] == 2 * d){

        colnames(res[[i]]) <- paste0(i, "_", rep( colnames(x), 2))

      } else {

        colnames(res[[i]]) <- i

      }
    }

    # Change objects to be coda compatible
    # Only apply to posterior samples
    if (attr(res[[i]], "type") == "sample"){

      # Differentiate between TVP and non TVP
      if (is.na(dim(res[[i]])[3]) == FALSE){

        # Create a sub list containing an mcmc object for each parameter in TVP case
        dat <- res[[i]]
        res[[i]] <- list()
        for (j in 1:dim(dat)[2]){
          res[[i]][[j]] <- as.mcmc(t(dat[, j, ]), start = niter - nburn, end = niter, thin = nthin)
          colnames(res[[i]][[j]]) <- paste0(i, "_", j, "_", 1:ncol(res[[i]][[j]]))

          # make it of class mcmc.tvp for custom plotting function
          class(res[[i]][[j]]) <- c("mcmc.tvp", "mcmc")

          attr(res[[i]][[j]], "type") <- "sample"

          # Imbue each mcmc.tvp object with index
          attr(res[[i]][[j]], "index") <- index
        }

        if (length(res[[i]]) == 1){
          res[[i]] <- res[[i]][[j]]
          attr(res[[i]][[j]], "index") <- index
        }

        # Make it of type 'sample' again
        attr(res[[i]], "type") <- "sample"

        # Rename
        if (dim(dat)[2] > 1){
          names(res[[i]]) <- colnames(dat)
        }


      } else {

        res[[i]] <- as.mcmc(res[[i]], start = niter - nburn, end = niter, thin = nthin)

      }
    }

    # Create summary of posterior
    if (is.list(res[[i]]) == FALSE & attr(res[[i]], "type") == "sample") {
      if (i != "theta_sr" & !(i == "sigma2" & sv == TRUE) & i != "beta") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){

          obj <- as.mcmc(x, start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
                          error = function(err) {
                            warning("Calculation of effective sample size failed for one or more variable(s). This can happen if the prior placed on the model induces extreme shrinkage.")
                            return(NA)
                          }, silent = TRUE)

          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = round(ESS)))
        }))
      } else if (i == "theta_sr") {
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){

          obj <- as.mcmc(abs(x), start = niter - nburn, end = niter, thin = nthin)
          ESS <- tryCatch(coda::effectiveSize(obj),
                          error = function(err) {
                            warning("Calculation of effective sample size failed for one or more variable(s). This can happen if the prior placed on the model induces extreme shrinkage.")
                            return(NA)
                          }, silent = TRUE)

          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = round(ESS)))
        }))
      }
    }
  }


  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }

  # add some attributes for the methods and plotting
  attr(res, "class") <- c("shrinkDTVP", "shrinkTVP")
  attr(res, "learn_a_xi") <- learn_a_xi
  attr(res, "learn_a_tau") <- learn_a_tau
  attr(res, "learn_c_xi") <- learn_c_xi
  attr(res, "learn_c_tau") <- learn_c_tau
  attr(res, "learn_kappa2_B") <- learn_kappa2_B
  attr(res, "learn_lambda2_B") <- learn_lambda2_B
  attr(res, "a_eq_c_xi") <- a_eq_c_xi
  attr(res, "a_eq_c_tau") <- a_eq_c_tau
  attr(res, "niter") <- niter
  attr(res, "nburn") <- nburn
  attr(res, "nthin") <- nthin
  attr(res, "sv") <- sv
  attr(res, "colnames") <-  colnames(x)
  attr(res, "index") <- index
  attr(res, "p") <- p
  attr(res, "mod_type") <- mod_type
  attr(res, "a_psi") <- a_psi
  attr(res, "c_psi") <- c_psi
  attr(res, "iid") <- iid
  attr(res, "shrink_inter") <- shrink_inter
  attr(res, "inter_column") <- inter_column




  return(res)
}
