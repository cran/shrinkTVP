#' Evaluate the one step ahead predictive density of a fitted TVP model
#'
#' \code{eval_pred_dens} evaluates the one-step ahead predictive density of a fitted TVP model resulting from a call to
#' shrinkTVP at the points supplied in \code{x}. For details on the approximation of the one-step ahead predictive density used,
#' see the vignette.
#'
#' @param x a real number or a vector of real numbers, taken to be the points at which the predictive density will be evaluated.
#' @param mod an object of class \code{shrinkTVP}, containing the fitted model for which the predictive density should be evaluated.
#' @param data_test a data frame with one row, containing the one step ahead covariates. The names of the covariates have to match the
#' names of the covariates used during model estimation in the call to \code{shrinkTVP}.
#' @param log a single logical value detrmining whether the density should be evaluated on the log scale or not.
#'
#' @return The value returned is a vector of length \code{length(x)}, containing the values of the predictive density evaluated
#' at the points supplied in \code{x}.
#'
#' @examples
#' \donttest{
#' # Simulate data
#'set.seed(123)
#'sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#'data <- sim$data
#'
#'# Estimate model
#'res <- shrinkTVP(y ~ x1 + x2, data = data[1:199, ])
#'
#'# Create sequence of x values where the density is to be evaluated
#'x_vals <- seq(0, 12, by = 0.1)
#'
#'# Evaluate density and plot
#'dens <- eval_pred_dens(x_vals, res, data[200, ])
#'plot(x_vals, dens, type = "l")
#'
#'# Add vertical line where true value of the one step ahead y lies
#'abline(v = data$y[200])
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#'@export
eval_pred_dens <- function(x, mod, data_test, log = FALSE){

  if (class(mod) != "shrinkTVP"){
    stop("mod has to be of class 'shrinkTVP'")
  }

  if (bool_input_bad(log)){
    stop("log has to be a single logical value")
  }

  if (is.data.frame(data_test) == FALSE || nrow(data_test) > 1){
    stop("data_test has to be a data frame with one row")
  }

  if (is.vector(x) == FALSE || any(sapply(x, numeric_input_bad_))){
    stop("x has to be a vector of numbers")
  }

  mf <- match.call(expand.dots = FALSE)
  mf$data <- data_test
  m <- match(x = "data", table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$formula <- mod$model$formula
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)

  # Remove the response from the data and formula
  mf$data <- mf$data[names(mf$data) != as.character(mf$formula)[2]]
  mf$formula[2] <- NULL

  mf <- eval(expr = mf, envir = parent.frame())
  mt <- attr(x = mf, which = "terms")
  # Create Matrix X with dummies and transformations
  x_test <- model.matrix(object = mt, data = mf)

  if (any(is.na(x_test))){
    stop("No NA values are allowed in covariates")
  }

  colnames(x_test)[colnames(x_test) == "(Intercept)"] <- "Intercept"

  if (attr(mod, "sv") == TRUE){
    sig2 <- mod$sigma2[, ncol(mod$sigma2)]
    sv_phi <- mod$sv_phi
    sv_mu <- mod$sv_mu
    sv_sigma2 <- mod$sigma2
  } else {
    sig2 <- mod$sigma2
    sv_phi = c(1)
    sv_mu = c(1)
    sv_sigma2 = c(1)
  }

  res <- pred_dens_mix_approx(theta_sr = as.matrix(mod$theta_sr),
                              beta_mean = as.matrix(mod$beta_mean),
                              sig2_samp = as.matrix(mod$sigma2),
                              sv = attr(mod, "sv"), sv_phi = sv_phi, sv_mu = sv_mu, sv_sigma2 = sv_sigma2,
                              x_test = x_test,
                              y_test = x,
                              chol_C_N_inv_samp = mod$LPDS_comp$chol_C_N_inv,
                              m_N_samp = mod$LPDS_comp$m_N,
                              M = nrow(mod$theta_sr), log = log)

  return(as.vector(res))

}

#' Calculate the Log Predictive Density Score for a fitted TVP model
#'
#'
#' \code{LPDS} calculates the one step ahead Log Predictive Density Score (LPDS) of a fitted TVP model resulting from a call to
#' shrinkTVP. For details on the approximation of the one-step ahead predictive density used, see the vignette.
#'
#' @param mod an object of class \code{shrinkTVP}, containing the fitted model for which the LPDS should be calculated.
#' @param data_test a data frame with one row, containing the one step ahead covariates and response. The names of the covariates
#' and the response have to match the names used during model estimation in the call to \code{shrinkTVP}.
#'
#' @return A real number equaling the calculated LPDS.
#'
#' @examples
#' \donttest{
#' # Simulate data
#'set.seed(123)
#'sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#'data <- sim$data
#'
#'# Estimate model
#'res <- shrinkTVP(y ~ x1 + x2, data = data[1:199, ])
#'
#'# Calculate LPDS
#' LPDS(res, data[200,])
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#'
#'@export
LPDS <- function(mod, data_test){

  if (class(mod) != "shrinkTVP"){
    stop("mod has to be of class 'shrinkTVP'")
  }

  if (is.data.frame(data_test) == FALSE || nrow(data_test) > 1){
    stop("data_test has to be a data frame with one row")
  }

  mf <- match.call(expand.dots = FALSE)
  mf$data <- data_test

  m <- match(x = "data", table = names(mf), nomatch = 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.pass
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(expr = mf, envir = parent.frame())
  # Create Vector y
  y <- model.response(mf, "numeric")

  return(as.numeric(eval_pred_dens(y, mod, data_test, log = TRUE)))
}

