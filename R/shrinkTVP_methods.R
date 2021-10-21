
#' Nicer printing of shrinkTVP objects
#'
#' @param x a \code{shrinkTVP} object.
#' @param ... Currently ignored.
#'
#' @return Called for its side effects and returns invisibly.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @export
print.shrinkTVP <- function(x, ...){
  ind <- attr(x, "index")
  cat(paste0("Object containing a fitted TVP model ", ifelse(attr(x, "sv"), "with stochastic volatility ", ""), "with:\n",
             " - ", formatC(length(attr(x, "colnames")), width = 7), " covariates", ifelse(attr(x, "p") > 0, paste0(", of which ", attr(x, "p"), " are AR terms"), ""), "\n",
             " - ", formatC(length(x$model$y), width = 7), " timepoints, running from ", min(ind), " to ", max(ind), "\n",
             " - ", formatC(attr(x, "niter"), width = 7), " MCMC draws\n",
             " - ", formatC(attr(x, "nburn"), width = 7), " burn-in\n",
             " - ", formatC(attr(x, "nthin"), width = 7), " thinning\n"))
  invisible(x)
}

#' @export
summary.shrinkTVP <- function(object, digits = 3, showprior = TRUE, ...) {

  # Check if digits is scalar, integer and positive
  if (is.scalar(digits) == FALSE |
      is.numeric(digits) == FALSE |
      digits %% 1 != 0 |
      digits < 0){
    stop("digits has to be a single, positive integer")
  }

  # Check if showprior is a logical
  if (bool_input_bad(showprior)){
    stop("showprior has to be a single logical value")
  }

  ret <- attributes(object)
  class(ret) <- c("summary.shrinkTVP")
  ret$priorvals <- object$priorvals
  ret$summaries <- object$summaries
  ret$types <- lapply(object, function(mod) return(attr(mod, "type")))
  ret$digits <- digits
  ret$showprior <- showprior
  ret
}

#' @method print summary.shrinkTVP
#' @export
print.summary.shrinkTVP <- function(x, ...) {
  cat("\nSummary of ", (x$niter - x$nburn), " MCMC draws after burn-in of ", x$nburn, ".\n", sep = "")

  if(x$showprior == TRUE){
    mod_type <- x$mod_type
    cat("\nPrior distributions:\n\n")

    if (mod_type == "double") {
      dist_parser("beta_mean | tau2 ~ Normal ( 0 , tau2 )", x)
      dist_parser("tau2 | a_tau , lambda2_B ~ Gamma ( a_tau , a_tau * lambda2_B / 2 )", x)
      dist_parser("a_tau | alpha_a_tau , beta_a_tau ~ Gamma ( alpha_a_tau , alpha_a_tau * beta_a_tau )", x)
      dist_parser("lambda2_B | e1 , e2 ~ Gamma ( e1 , e2 )", x)
      cat("\n")
      dist_parser("theta_sr | xi2 ~ Normal ( 0 , xi2 )", x)
      dist_parser("xi2 | a_xi , kappa2_B ~ Gamma ( a_xi , a_xi * kappa2_B / 2 )", x)
      dist_parser("a_xi | alpha_a_xi , beta_a_xi ~ Gamma ( alpha_a_xi , alpha_a_xi * beta_a_xi )", x)
      dist_parser("kappa2_B | d1 , d2 ~ Gamma ( d1 , d2 )", x)
    } else  if (mod_type == "triple") {
      if (x$a_eq_c_tau == FALSE) {
        dist_parser("beta_mean | a_tau , c_tau , tau2 , lambda2 , lambda2_B ~ Normal ( 0 , 2 *  c_tau * tau2 / ( a_tau * lambda2 * lambda2_B ) )", x)
        dist_parser("tau2 | a_tau  ~ Gamma ( a_tau , 1)", x)
        dist_parser("lambda2 | c_tau ~ Gamma ( c_tau , 1)", x)
        dist_parser("lambda2_B / 2 | a_tau , c_tau ~ F ( 2 * a_tau , 2 * c_tau )", x)
        dist_parser("2 * a_tau | alpha_a_tau , beta_a_tau ~ Beta ( alpha_a_tau , beta_a_tau )", x)
        dist_parser("2 * c_tau | alpha_c_tau , beta_c_tau ~ Beta ( alpha_c_tau , beta_c_tau )", x)
      } else {
        dist_parser("beta_mean | a_tau , tau2 , lambda2 , lambda2_B ~ Normal ( 0 , 2 *  a_tau * tau2 / ( a_tau * lambda2 * lambda2_B ) )", x)
        dist_parser("tau2 | a_tau  ~ Gamma ( a_tau , 1)", x)
        dist_parser("lambda2 | a_tau  ~ Gamma ( a_tau , 1 )", x)
        dist_parser("lambda2_B / 2 | a_tau ~ F ( 2 * a_tau , 2 * a_tau )", x)
        dist_parser("2 * a_tau = 2 * c_tau | alpha_a_tau , beta_a_tau ~ Beta ( alpha_a_tau , beta_a_tau )", x)
      }

      cat("\n")
      if (x$a_eq_c_xi == FALSE) {
        dist_parser("theta_sr | a_xi , c_xi , xi2 , kappa2 , kappa2_B ~  Normal ( 0 , 2 *  c_xi * xi2 / ( a_xi * kappa2 * kappa2_B ) )", x)
        dist_parser("xi2 | a_xi ~ Gamma ( a_xi , 1 )", x)
        dist_parser("kappa2 | c_xi ~ Gamma ( c_xi , 1 )", x)
        dist_parser("kappa2_B | a_xi , c_xi ~ F ( 2 * a_xi , 2 * c_xi )", x)
        dist_parser("2 * a_xi | alpha_a_xi , beta_a_xi ~ Beta ( alpha_a_xi , beta_a_xi )", x)
        dist_parser("2 * c_xi | alpha_c_xi , beta_c_xi ~ Beta ( alpha_c_xi , beta_c_xi )", x)
      } else {
        dist_parser("theta_sr | a_xi , xi2 , kappa2 , kappa2_B ~  Normal ( 0 , 2 *  a_xi * xi2 / ( a_xi * kappa2 * kappa2_B ) )", x)
        dist_parser("xi2 | a_xi ~ Gamma ( a_xi , 1 )", x)
        dist_parser("lambda2 | a_tau ~ Gamma ( a_tau , 1 )", x)
        dist_parser("kappa2_B | a_xi ~ F ( 2 * a_xi , 2 * a_xi )", x)
        dist_parser("2 * a_xi = 2 * c_xi | alpha_a_xi , beta_a_xi ~ Beta ( alpha_a_xi , beta_a_xi )", x)
      }
    } else if (mod_type == "ridge") {
      dist_parser("beta_mean | lambda2_B ~ Normal ( 0 , 2 / lambda2_B )", x)
      dist_parser("theta_sr | kappa2_B ~ Normal ( 0 , 2 / kappa2_B )", x)
    }

    # Prior distributions for sigma2
    if (x$sv == TRUE){
      cat("\nPrior distributions on stochastic volatility parameters:\n")
      dist_parser("sv_mu ~ Normal ( bmu , Bmu )", x)
      dist_parser("( sv_phi + 1 ) / 2 ~ Beta ( a0_sv , b0_sv )", x)
      dist_parser("sv_sigma2 ~ Gamma ( 1/2 , 1 / ( 2 * Bsigma_sv ) )", x)
    } else {
      cat("\n")
      dist_parser(" sigma2 | C0 ~ GammaInv ( c0 , C0 )", x)
      dist_parser(" C0 ~ Gamma ( g0 , G0 )", x)
    }

    cat("\n")

  }

  # Posterior summaries
  cat("\nStatistics of posterior draws of parameters (thinning = ", x$nthin, "):\n\n", sep = "")

  # The posterior summaries are printed by printing a dataframe where all NA values are
  # set to "", thereby not showing up in the console. This is done to make sure everything
  # is nicely aligned
  ind <- 1
  for (i in x$summaries){
    if (ind == 1){
      post_sum <- rbind(round(i, x$digits), matrix(NA, ncol = ncol(i), nrow = 1))
    } else {
      post_sum <- rbind(post_sum, round(i, x$digits), matrix(NA, ncol = ncol(i), nrow = 1))
    }

    ind <- ind + 1
  }
  post_sum <- as.data.frame(cbind(rownames(post_sum), post_sum), stringsAsFactors = FALSE, row.names = FALSE)
  post_sum[is.na(post_sum)] <- ""

  # Adjust column names
  colnames(post_sum) <- c("param", "mean", "sd", "median", "HPD 2.5%", "HPD 97.5%", "ESS")

  # The summaries of theta_sr are based on the absolute value, this has to be reflected in the param name
  post_sum$param <- ifelse(grepl("theta_sr", post_sum$param), paste0("abs(", post_sum$param, ")"), post_sum$param)

  # Finally print it
  print(post_sum, row.names = FALSE, right = FALSE)

  invisible(x)
}

#' Calculate fitted historical values for an estimated TVP model
#'
#' Calculates the fitted values for an estimated TVP model, i.e. \eqn{X_t'\beta_t}.
#' Note that in contrast to \code{\link{predict.shrinkTVP}} this does not include the error term.
#'
#' @param object A \code{shrinkTVP} object
#' @param ... Currently ignored.
#'
#' @return An object of class \code{shrinkTVP_fitted}
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' \donttest{
#'
#' # Generate synthetic data
#' sim <- simTVP()
#'
#' # Estimate a model
#' res <- shrinkTVP(y ~ x1 + x2, sim$data)
#'
#' # Calculate fitted values
#' fitted <- fitted(res)
#'
#' # Visualize
#' plot(fitted)
#' lines(sim$data$y, col = "forestgreen")
#' }
#' @family prediction functions
#' @export
fitted.shrinkTVP <- function(object, ...){

  fitted <- calc_fitted(object$model$y, object$model$x, object$beta)
  class(fitted) <- c("shrinkTVP_fitted", "mcmc.tvp")
  attr(fitted, "index") <- attr(object, "index")

  return(fitted)
}

#' Calculate residuals for an estimated TVP model
#'
#' Calculates the residuals for an estimated TVP model, i.e. \eqn{y_t - X_t'\beta_t}.
#'
#' @param object a \code{shrinkTVP} object.
#' @param ... Currently ignored.
#'
#' @return An object of class \code{shrinkTVP_resid}
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' \donttest{
#'
#' # Generate synthetic data
#' sim <- simTVP(N = 300)
#'
#' # Estimate a model
#' res <- shrinkTVP(y ~ x1 + x2, sim$data)
#'
#' # Calculate residuals
#' resids <- residuals(res)
#'
#' # Visualize
#' plot(resids)
#' }
#' @family prediction functions
#' @export
residuals.shrinkTVP <- function(object, ...){

  fitted <- fitted(object)
  resids <- t(- t(fitted) + c(object$model$y))
  class(resids) <- c("shrinkTVP_resid", "mcmc.tvp")
  attr(resids, "index") <- attr(object, "index")

  return(resids)
}

#' Calculate predicted historical values for an estimated TVP model
#'
#' Calculates the predicted past values for an estimated TVP model, i.e. \eqn{X_t'\beta_t + \epsilon_t}.
#' Note that in contrast to \code{\link{fitted.shrinkTVP}} this includes the error term.
#'
#' @param object a \code{shrinkTVP} object
#' @param ... Currently ignored.
#'
#' @return An object of class \code{shrinkTVP_pred}.
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @examples
#' \donttest{
#'
#' # Generate synthetic data
#' sim <- simTVP(N = 300)
#'
#' # Estimate a model
#' res <- shrinkTVP(y ~ x1 + x2, sim$data)
#'
#' # Calculate predicted values
#' pred <- predict(res)
#'
#' # Visualize
#' plot(pred)
#' lines(sim$data$y, col = "forestgreen")
#' }
#' @family prediction functions
#' @export
predict.shrinkTVP <- function(object, ...){

  nsave <- nrow(object$sigma2)
  nT <- length(attr(object, "index"))

  fitted <- calc_fitted(object$model$y, object$model$x, object$beta)
  pred <- fitted + matrix(rnorm(nT * nsave, 0, sqrt(object$sigma2)), nsave, nT, byrow = FALSE)
  class(pred) <- c("shrinkTVP_pred", "mcmc.tvp")
  attr(pred, "index") <- attr(object, "index")

  return(pred)
}
