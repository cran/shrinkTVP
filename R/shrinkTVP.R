
#' Markov Chain Monte Carlo (MCMC) for time-varying parameter models with shrinkage
#'
#' \code{shrinkTVP} samples from the joint posterior distribution of the parameters of a time-varying
#' parameter model with shrinkage, potentially including stochastiv volatility (SV), and returns the MCMC draws.
#'
#' For details concerning the algorithm please refer to the paper by Bitto and Frühwirth-Schnatter (2019).
#'
#' @param formula an object of class "formula": a symbolic representation of the model, as in the
#' function \code{lm}. For details, see \code{\link{formula}}.
#' @param data \emph{optional} data frame containing the response variable and the covariates. If not found in \code{data},
#' the variables are taken from \code{environment(formula)}, typically the environment from which \code{shrinkTVP}
#' is called. No \code{NA}s are allowed in the response variable and the covariates.
#' @param niter positive integer, indicating the number of MCMC iterations
#' to perform, including the burn-in. Has to be larger than or equal to \code{nburn} + 2. The default value is 10000.
#' @param nburn non-negative integer, indicating the number of iterations discarded
#' as burn-in. Has to be smaller than or equal to \code{niter} - 2. The default value is \code{round(niter / 2)}.
#' @param nthin positive integer, indicating the degree of thinning to be performed. Every \code{nthin} draw is kept and returned.
#' The default value is 1, implying that every draw is kept.
#' @param learn_a_xi logical value indicating whether to learn a_xi, the local adaptation parameter of the state variances.
#' The default value is \code{TRUE}.
#' @param learn_a_tau logical value indicating whether to learn a_tau, the local adaptation parameter of the mean of the
#' initial values of the states. The default value is \code{TRUE}.
#' @param a_xi positive, real number, indicating the (fixed) value for a_xi. Ignored if
#' \code{learn_a_xi} is \code{TRUE}. The default value is 0.1.
#' @param a_tau positive, real number, indicating the (fixed) value for a_tau. Ignored if
#' \code{learn_a_tau} is \code{TRUE}. The default value is 0.1.
#' @param learn_kappa2 logical value indicating whether to learn kappa2, the global level of shrinkage for
#' the state variances. The default value is \code{TRUE}.
#' @param learn_lambda2 logical value indicating whether to learn the lambda^2 parameter,
#' the global level of shrinkage for the mean of the initial values of the states. The default value is \code{TRUE}.
#' @param kappa2 positive, real number, indicating the (fixed) value for kappa2. Ignored if
#' \code{learn_kappa2} is \code{TRUE}. The default value is 20.
#' @param lambda2 positive, real number, indicating the (fixed) value for lambda2. Ignored if
#' \code{learn_lambda2} is \code{TRUE}. The default value is 20.
#' @param hyperprior_param \emph{optional} named list containing hyperparameter values.
#' Not all have to be supplied, with those missing being replaced by the default values.
#' Any list elements that are misnamed will be ignored and a warning will be thrown.
#' All hyperparameter values have to be positive, real numbers. The following hyperparameters can be
#' supplied:
#' \itemize{
#' \item \code{c0}: The default value is 2.5.
#' \item \code{g0}: The default value is 5.
#' \item \code{G0}: The default value is 5 / (2.5 - 1).
#' \item \code{cp}: The default value is 1.
#' \item \code{np}: The default value is 20.
#' \item \code{e1}: The default value is 0.001.
#' \item \code{e2}: The default value is 0.001.
#' \item \code{d1}: The default value is 0.001.
#' \item \code{d2}: The default value is 0.001.
#' \item \code{b_xi}: The default value is 10.
#' \item \code{b_tau}: The default value is 10.
#' \item \code{nu_xi}: The default value is 5.
#' \item \code{nu_tau}: The default value is 5.
#' }
#' @param c_tuning_par_xi positive, real number. Determines the standard deviation of the proposal
#' distribution for the Metropolis Hastings step for a_xi. Ignored if \code{learn_a_xi} is \code{FALSE}. The default
#' value is 1.
#' @param c_tuning_par_tau positive, real number. Determines the standard deviation of the proposal
#' distribution for the Metropolis Hastings step for a_tau. Ignored if \code{learn_a_tau} is \code{FALSE}. The default
#' value is 1.
#' @param display_progress logical value indicating whether the progress bar and other informative output should be
#' displayed. The default value is \code{TRUE}.
#' @param ret_beta_nc logical value indicating whether to output the non-centered states in addition to the centered ones.
#' The default value is \code{FALSE}.
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
#' @param LPDS logical value indicating whether the one step ahead log predictive density score should be returned.
#' If \code{LPDS} is \code{TRUE}, both \code{x_test} and \code{y_test} have to be supplied. The default value is \code{FALSE}.
#' @param y_test \emph{optional} real number. If \code{LPDS} is \code{TRUE} this value will be taken to be the true one step
#' ahead value and will be used to calculate the log predictive density score. Ignored if \code{LPDS} is \code{FALSE}.
#' @param x_test \emph{optional} object that can be coerced to a 1 x d data frame containing the covariates for the
#' calculation of the one step ahead log predictive density score. The column names have to exactly match the
#' names of the covariates provided in the formula or an error will be thrown. An intercept will be added to
#' \code{x_test} if an intercept was added to the covariates by the formula interface. Ignored if \code{LPDS} is \code{FALSE}.
#'
#' @return The value returned is a list object of class \code{shrinkTVP_res} containing
#' \item{\code{sigma2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of \code{sigma2}.
#' If \code{sv} is \code{TRUE}, \code{sigma2} is additionally an \code{mcmc.tvp} object.}
#' \item{\code{theta_sr}}{an \code{mcmc} object containing the parameter draws from the posterior distribution of the square root of theta.}
#' \item{\code{beta_mean}}{an \code{mcmc} object containing the parameter draws from the posterior distribution of beta_mean.}
#' \item{\code{beta_nc}}{\emph{(optional)} \code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior
#' distribution of the non-centered states, one for each covariate. In the case that there is only one covariate, this becomes just
#' a single \code{mcmc.tvp} object. Not returned if \code{ret_beta_nc} is \code{FALSE}.}
#' \item{\code{beta}}{\code{list} object containing an \code{mcmc.tvp} object for the parameter draws from the posterior distribution of the centered
#' states, one for each covariate. In the case that there is only one covariate, this becomes just a single \code{mcmc.tvp} object.}
#' \item{\code{xi2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of xi2.}
#' \item{\code{a_xi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_xi.
#' Not returned if \code{learn_a_xi} is \code{FALSE}.}
#' \item{\code{a_xi_acceptance}}{\emph{(optional)} \code{list} object containing acceptance statistics for the Metropolis Hastings algorithm for
#' a_xi. Not returned if \code{learn_a_xi} is \code{FALSE}.}
#' \item{\code{tau2}}{\code{mcmc} object containing the parameter draws from the posterior distribution of tau2.}
#' \item{\code{a_tau}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of a_tau.
#' Not returned if \code{learn_a_tau} is \code{FALSE}.}
#' \item{\code{a_tau_acceptance}}{\emph{(optional)} \code{list} containing acceptance statistics for the Metropolis Hastings algorithm for
#' a_tau. Not returned if \code{learn_a_tau} is \code{FALSE}.}
#' \item{\code{kappa2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of kappa2.
#' Not returned if \code{learn_kappa2} is \code{FALSE}.}
#' \item{\code{lambda2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of lambda2.
#' Not returned if \code{learn_lambda2} is \code{FALSE}.}
#' \item{\code{C0}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of C0.
#' Not returned if \code{sv} is \code{TRUE}.}
#' \item{\code{sv_mu}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the mu
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{sv_phi}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the phi
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{sv_sigma2}}{\emph{(optional)} \code{mcmc} object containing the parameter draws from the posterior distribution of the sigma2
#' parameter for the stochastic volatility model on the errors. Not returned if \code{sv} is \code{FALSE}.}
#' \item{\code{priorvals}}{\code{list} object containing hyperparameter values of the prior distributions, as specified by the user.}
#' \item{\code{model}}{\code{list} object containing the model matrix and model response used.}
#' \item{\code{summaries}}{\code{list} object containing a collection of summary statistics of the posterior draws.}
#' \item{\code{LPDS}}{\emph{(optional)} value of the log predictive density score, calculated with \code{y_test} and \code{x_test}.}
#' To display the output, use \code{plot} and \code{summary}. The \code{summary} method displays the specified prior values stored in
#' \code{priorvals} and the posterior summaries stored in \code{summaries}, while the \code{plot} method calls \code{coda}'s \code{plot.mcmc}
#' or the \code{plot.mcmc.tvp} method. Furthermore, all functions that can be appplied to \code{coda::mcmc} objects
#' (e.g. \code{coda::acfplot}) can be applied to all output elements that are \code{coda} compatible.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @seealso \code{\link{plot.shrinkTVP_res}}, \code{\link{plot.mcmc.tvp}}
#' @references Bitto, A., & Frühwirth-Schnatter, S. (2019). "Achieving shrinkage in a time-varying parameter model framework."
#' \emph{Journal of Econometrics}, 210(1), 75-97. <doi:10.1016/j.jeconom.2018.11.006>
#' @examples
#' \donttest{
#'
#' ## Example 1, learn everything
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' res <- shrinkTVP(y ~ x1 + x2, data = data)
#' # summarize output
#' summary(res)
#'
#'
#' ## Example 2, hierarchical Bayesian Lasso
#' res <- shrinkTVP(y ~ x1 + x2, data = data,
#'                 learn_a_xi = FALSE, learn_a_tau = FALSE,
#'                 a_xi = 1, a_tau = 1)
#'
#'
#' ## Example 3, non-hierarchical Bayesian Lasso
#' res <- shrinkTVP(y ~ x1 + x2, data = data,
#'                 learn_a_xi = FALSE, learn_a_tau = FALSE,
#'                 a_xi = 1, a_tau = 1,
#'                 learn_kappa2 = FALSE, learn_lambda2 = FALSE)
#'
#'
#' ## Example 4, adding stochastic volatility
#' res <- shrinkTVP(y ~ x1 + x2, data = data,
#'                 sv = TRUE)
#'
#'
#' ## Example 4, changing some of the default hyperparameters
#' res <- shrinkTVP(y ~ x1 + x2, data = data,
#'                 hyperprior_param = list(cp = 5,
#'                                         nu_xi = 10))
#' }
#'
#' @export
shrinkTVP <- function(formula,
                      data,
                      niter = 10000,
                      nburn = round(niter / 2),
                      nthin = 1,
                      learn_a_xi = TRUE,
                      learn_a_tau = TRUE,
                      a_xi = 0.1,
                      a_tau = 0.1,
                      learn_kappa2 = TRUE,
                      learn_lambda2 = TRUE,
                      kappa2 = 20,
                      lambda2 = 20,
                      hyperprior_param,
                      c_tuning_par_xi = 1,
                      c_tuning_par_tau = 1,
                      display_progress = TRUE,
                      ret_beta_nc = FALSE,
                      sv = FALSE,
                      sv_param,
                      LPDS  = FALSE,
                      y_test,
                      x_test){


  # Input checking ----------------------------------------------------------


  # default hyperparameter values
  default_hyper <- list(c0 = 2.5,
                         g0 = 5,
                         G0 = 5 / (2.5 - 1),
                         cp = 1,
                         np = 20,
                         e1 = 0.001,
                         e2 = 0.001,
                         d1 = 0.001,
                         d2 = 0.001,
                         b_xi = 10,
                         b_tau = 10,
                         nu_xi = 5,
                         nu_tau = 5)

  # default sv params
  default_hyper_sv <- list(Bsigma_sv = 1,
                            a0_sv = 5,
                            b0_sv = 1.5,
                            bmu = 0,
                            Bmu = 1)

  # Change hyperprior values if user overwrites them
  if (missing(hyperprior_param)){
    hyperprior_param <- default_hyper
  } else {

    # Check that hyperprior_param and sv_param are a list
    if (is.list(hyperprior_param) == FALSE | is.data.frame(hyperprior_param)){
      stop("hyperprior_param has to be a list")
    }

    stand_nam <- names(default_hyper)
    user_nam <- names(hyperprior_param)

    # Give out warning if an element of the parameter list is misnamed
    if (any(!user_nam %in% stand_nam)){
      wrong_nam <- user_nam[!user_nam %in% stand_nam]
      warning(paste0(paste(wrong_nam, collapse = ", "),
                     ifelse(length(wrong_nam) == 1, " has", " have"),
                     " been incorrectly named in hyperprior_param and will be ignored"),
              immediate. = TRUE)
    }

    # Merge users' and default values and ignore all misnamed values
    missing_param <- stand_nam[!stand_nam %in% user_nam]
    hyperprior_param[missing_param] <- default_hyper[missing_param]
    hyperprior_param <- hyperprior_param[stand_nam]
  }


  # Same procedure for sv_param
  if (missing(sv_param) | sv == FALSE){
    sv_param <- default_hyper_sv
  } else {

    if (is.list(sv_param) == FALSE | is.data.frame(sv_param)){
      stop("sv_param has to be a list")
    }

    stand_sv_nam <- names(default_hyper_sv)
    user_sv_nam <- names(sv_param)

    # Give out warning if an element of the parameter list is misnamed
    if (any(!user_sv_nam %in% stand_sv_nam)){
      wrong_nam <- user_sv_nam[!user_sv_nam %in% stand_sv_nam]
      warning(paste0(paste(wrong_nam, collapse = ", "),
                     ifelse(length(wrong_nam) == 1, " has", " have"),
                     " been incorrectly named in sv_param and will be ignored"),
              immediate. = TRUE)
    }

    # Merge users' and default values and ignore all misnamed values
    missing_sv_param <- stand_sv_nam[!stand_sv_nam %in% user_sv_nam]
    sv_param[missing_sv_param] <- default_hyper_sv[missing_sv_param]
    sv_param <- sv_param[stand_sv_nam]
  }

  # Check if all numeric inputs are correct
  to_test_num <- list(lambda2 = lambda2,
                   kappa2 = kappa2,
                   a_xi = a_xi,
                   a_tau = a_tau,
                   c_tuning_par_xi = c_tuning_par_xi,
                   c_tuning_par_tau = c_tuning_par_tau)

  if (missing(hyperprior_param) == FALSE){
    to_test_num <- c(to_test_num, hyperprior_param)
  }

  if (missing(sv_param) == FALSE){
    to_test_num <- c(to_test_num, sv_param[names(sv_param) != "bmu"])
  }

  bad_inp <- sapply(to_test_num, numeric_input_bad)


  if (any(bad_inp)){
    stand_names <- c(names(default_hyper), names(default_hyper_sv), "lambda2", "kappa2", "a_xi", "a_tau", "c_tuning_par_xi", "c_tuning_par_tau")
    bad_inp_names <- names(to_test_num)[bad_inp]
    bad_inp_names <- bad_inp_names[bad_inp_names %in% stand_names]
    stop(paste0(paste(bad_inp_names, collapse = ", "),
                ifelse(length(bad_inp_names) == 1, " has", " have"),
                " to be a single, positive number"))
  }

  # Check bmu seperately
  if (!is.numeric(sv_param$bmu) | !is.scalar(sv_param$bmu)){
    stop("bmu has to be a single number")
  }

  # Check if all integer inputs are correct
  to_test_int <- list(niter = niter,
                   nburn = nburn,
                   nthin = nthin)
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
  to_test_bool <- list(learn_lambda2 = learn_lambda2,
                       learn_kappa2 = learn_kappa2,
                       learn_a_xi = learn_a_xi,
                       learn_a_tau = learn_a_tau,
                       display_progress = display_progress,
                       sv = sv,
                       LPDS  = LPDS,
                       ret_beta_nc = ret_beta_nc)

  bad_bool_inp <- sapply(to_test_bool, bool_input_bad)

  if (any(bad_bool_inp)){
    bad_inp_names <- paste(names(to_test_bool)[bad_bool_inp], collapse = ", ")
    stop(paste0(bad_inp_names,
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

  colnames(x)[colnames(x) == "(Intercept)"] <- "Intercept"

  d <- dim(x)[2]
  a0 <- rep(0, 2 * d)
  store_burn <- FALSE

  # Check y_test and x_test

  if (LPDS == TRUE){
    if (missing(x_test)){
      stop("x_test is missing")
    }
    if (missing(y_test)){
      stop("y_test is missing")
    }

  }

  if (missing(y_test)){
    y_test <- NA
  } else {
    if (LPDS == TRUE) {
      if(is.scalar(y_test) == FALSE | is.numeric(y_test) == FALSE){
        stop("y_test has to be a single number")
      }
    } else {
      y_test <- NA
    }
  }

  if (missing(x_test)){
    x_test <- NA
  } else {
    if (LPDS == TRUE){


      x_test <- tryCatch(as.data.frame(x_test),
                         error = function(e) stop("x_test could not be coerced to data frame"))

      # Check that x_test only has one row
      if (nrow(x_test) > 1){
        stop("The coerced x_test data frame had more than one row")
      }

      # Check that all colnames (except 'Intercept') are present in x_test
      missing <- colnames(x)[!colnames(x) %in% colnames(x_test) & colnames(x) != "Intercept"]
      if (length(missing) > 0){
        missing_names <- paste(missing, collapse = ", ")
        stop(paste0(missing_names,
                    ifelse(length(missing_names) == 1, " has", " have"),
                    " to be present in x_test"))
      }

      # Throw warning if user misnamed columns
      misnamed <- colnames(x_test)[!colnames(x_test) %in% colnames(x)]
      if (length(misnamed) > 0){
        misnamed_x_test <- paste(misnamed, collapse = ", ")
        warning(paste0(misnamed_x_test,
                       ifelse(length(misnamed) > 1, " have", " has"),
                       " been incorrectly named in x_test and will be ignored"),
                immediate. = TRUE)
      }

      # Bring variables into correct order
      x_test <- x_test[,colnames(x)[colnames(x) != "Intercept"]]

      # Turn into matrix, as this is what the C++ code expects
      x_test <- as.matrix(x_test)

      # Add Intercept if the formula added one
      if ("Intercept" %in% colnames(x)) {
        x_test <- cbind(1, x_test)
      }

    } else {
      x_test <- NA
    }
  }


  # Run sampler -------------------------------------------------------------


  runtime <- system.time({
    suppressWarnings({
      res <- do_shrinkTVP(y,
                          x,
                          a0,
                          niter,
                          nburn,
                          nthin,
                          hyperprior_param$c0,
                          hyperprior_param$g0,
                          hyperprior_param$G0,
                          hyperprior_param$cp,
                          hyperprior_param$np,
                          hyperprior_param$d1,
                          hyperprior_param$d2,
                          hyperprior_param$e1,
                          hyperprior_param$e2,
                          learn_lambda2,
                          learn_kappa2,
                          lambda2,
                          kappa2,
                          learn_a_xi,
                          learn_a_tau,
                          a_xi,
                          a_tau,
                          c_tuning_par_xi,
                          c_tuning_par_tau,
                          hyperprior_param$b_xi,
                          hyperprior_param$b_tau,
                          hyperprior_param$nu_xi,
                          hyperprior_param$nu_tau,
                          display_progress,
                          ret_beta_nc,
                          store_burn,
                          sv,
                          sv_param$Bsigma_sv,
                          sv_param$a0_sv,
                          sv_param$b0_sv,
                          sv_param$bmu,
                          sv_param$Bmu,
                          y_test,
                          x_test,
                          LPDS)
    })
  })

  # Throw an error if the sampler failed
  if (res$success_vals$success == FALSE){
    stop(paste0("The sampler failed at iteration ",
                res$success_vals$fail_iter,
                " while trying to ",
                res$success_vals$fail, ". ",
                "Try rerunning the model. ",
                "If the sampler fails again, try changing the prior to be more informative. ",
                "If the problem still persists, please contact the maintainer: ",
                maintainer("shrinkTVP")))
  } else {
    res$success_vals <- NULL
  }


  # Post process sampler results --------------------------------------------


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
  if (ret_beta_nc == FALSE){
    res[["beta_nc"]] <- NULL
  }

  if (sv == TRUE){
    res[["C0"]] <- NULL
  } else {
    res[["sv_mu"]] <- NULL
    res[["sv_phi"]] <- NULL
    res[["sv_sigma2"]] <- NULL
  }

  if (learn_kappa2 == FALSE){
    res[["kappa2"]] <- NULL
  }

  if (learn_lambda2 == FALSE){
    res[["lambda2"]] <- NULL
  }

  if (learn_a_xi == FALSE){
    res[["a_xi"]] <- NULL
    res[["a_xi_acceptance"]] <- NULL
  }

  if (learn_a_tau == FALSE){
    res[["a_tau"]] <- NULL
    res[["a_tau_acceptance"]] <- NULL
  }

  if (LPDS == FALSE){
    res[["LPDS"]] <- NULL
  }

  # Create object to hold prior values
  priorvals <- c()

  priorvals["np"] <- hyperprior_param$np
  priorvals["cp"] <- hyperprior_param$cp

  if (learn_a_tau == TRUE){
    priorvals["b_tau"] <- hyperprior_param$b_tau
    priorvals["nu_tau"] <- hyperprior_param$nu_tau
  } else {
    priorvals["a_tau"] <- a_tau
  }

  if (learn_a_xi == TRUE){
    priorvals["b_xi"] <- hyperprior_param$b_xi
    priorvals["nu_xi"] <- hyperprior_param$nu_xi
  } else {
    priorvals["a_xi"] <- a_xi
  }


  if (learn_lambda2 == TRUE){
    priorvals["e1"] <- hyperprior_param$e1
    priorvals["e2"] <- hyperprior_param$e2
  } else {
    priorvals["lambda2"] <- lambda2
  }

  if (learn_kappa2 == TRUE){
    priorvals["d1"] <- hyperprior_param$d1
    priorvals["d2"] <- hyperprior_param$d2
  } else {
    priorvals["kappa2"] <- kappa2
  }

  if (sv == TRUE){
    priorvals["Bsigma_sv"] <- sv_param$Bsigma_sv
    priorvals["a0_sv"] <- sv_param$a0_sv
    priorvals["b0_sv"] <- sv_param$b0_sv
    priorvals["bmu"] <- sv_param$bmu
    priorvals["Bmu"] <- sv_param$Bmu
  } else {
    priorvals["g0"] <- hyperprior_param$g0
    priorvals["G0"] <- hyperprior_param$G0
    priorvals["c0"] <- hyperprior_param$c0
  }

  res$priorvals <- priorvals

  # Add data to output
  res[["model"]] <- list()
  res$model$x <- x
  res$model$y <- y

  res$summaries <- list()

  # add attributes to the individual if they are distributions or individual statistics
  nsave <- ifelse(store_burn, floor(niter/nthin), floor((niter - nburn)/nthin))
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
        }

        if (length(res[[i]]) == 1){
          res[[i]] <- res[[i]][[j]]
        }

        # Make it of type 'sample' again
        attr(res[[i]], "type") <- "sample"

        # Rename
        names(res[[i]]) <- colnames(dat)


      } else {

        res[[i]] <- as.mcmc(res[[i]], start = niter - nburn, end = niter, thin = nthin)

      }
    }

    # Create summary of posterior
    if (is.list(res[[i]]) == FALSE & attr(res[[i]], "type") == "sample"){
      if (i != "theta_sr" & !(i == "sigma2" & sv == TRUE)){
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){
          obj <- as.mcmc(x, start = niter - nburn, end = niter, thin = nthin)
          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = effectiveSize(obj)))
        }))
      } else if (i == "theta_sr"){
        res$summaries[[i]] <- t(apply(res[[i]], 2, function(x){
          obj <- as.mcmc(abs(x), start = niter - nburn, end = niter, thin = nthin)
          return(c("mean" = mean(obj),
                   "sd" = sd(obj),
                   "median" = median(obj),
                   "HPD" = HPDinterval(obj)[c(1, 2)],
                   "ESS" = effectiveSize(obj)))
        }))
      }
    }
  }


  if (display_progress == TRUE) {
    cat("Done!\n", file = stderr())
  }

  # add some attributes for the methods and plotting
  attr(res, "class") <- "shrinkTVP_res"
  attr(res, "learn_a_xi") <- learn_a_xi
  attr(res, "learn_a_tau") <- learn_a_tau
  attr(res, "learn_kappa2") <- learn_kappa2
  attr(res, "learn_lambda2") <- learn_lambda2
  attr(res, "niter") <- niter
  attr(res, "nburn") <- nburn
  attr(res, "nthin") <- nthin
  attr(res, "sv") <- sv
  attr(res, "colnames") <-  colnames(x)



  return(res)
}
