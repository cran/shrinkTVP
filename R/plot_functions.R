#' Graphical summary of posterior distribution for a time-varying parameter
#'
#' \code{plot.mcmc.tvp} plots empirical posterior quantiles for a time-varying parameter.
#'
#' @param x \code{mcmc.tvp} object
#' @param probs vector of boundaries for credible intervals to plot for each point in time, with values in [0,1].
#' The largest and smallest value form the outermost credible interval, the second smallest and second largest the second outermost and so forth.
#' The default value is \code{c(0.025, 0.25, 0.75, 0.975)}. Note that there have to be the same number of probs
#' < 0.5 as there are > 0.5.
#' @param shaded single logical value or a vector of logical values, indicating whether or not to shade the area between the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine if the area between the outermost credible interval
#' is shaded, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of quantile pairs. The default value is \code{TRUE}.
#' @param quantlines single logical value or a vector of logical values, indicating whether or not to draw borders along the pointwise
#' credible intervals. If a vector is given, the first value given is used to determine whether the outermost credible interval
#' is marked by lines, the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the
#' number of credible intervals. The default value is \code{FALSE}.
#' @param shadecol single character string or a vector of character strings. Determines the color of the shaded areas that represent
#' the credible intervals. If a vector is given, the first color given is used for the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{shaded = FALSE}. The default value is \code{"skyblue"}.
#' @param shadealpha real number between 0 and 1 or a vector of real numbers between 0 and 1.
#' Determines the level of transparency of the shaded areas that represent
#' the credible intervals. If a vector is given, the first value
#' given is used for the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion
#' if the vector is shorter than the number of shaded areas. Has no effect if \code{shaded = FALSE}.
#' The default value is \code{0.5}.
#' @param quantlty either a single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"} or a vector containing these. Determines the line type of the borders
#' drawn around the shaded areas that represent the credible intervals. Note that if a vector is supplied the elements have to either
#' be all integers or all character strings. If a vector is given, the first value given is used for the outermost area, the second for
#' the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number of shaded areas.
#' Has no effect if \code{quantlines = FALSE}. The default value is \code{2}.
#' @param quantcol single character string or a vector of character strings. Determines the color of the borders drawn around the shaded
#' areas that represent the credible intervals. If a vector is given, the first color given is used for borders of the outermost area,
#' the second for the second outermost and so forth. Recycled in the usual fashion if the vector is shorter than the number
#' of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{"black"}.
#' @param quantlwd single real, positive number or a vector of real, positive numbers. Determines the line width of the borders
#' drawn around the shaded areas that represent the credible intervals. If a vector is given, the first number given is used for
#' the borders of the outermost area, the second for the second outermost and so forth. Recycled in the usual fashion if the vector
#' is shorter than the number of shaded areas. Has no effect if \code{quantlines = FALSE}. The default value is \code{0.5}.
#' @param drawzero single logical value determining whether to draw a horizontal line at zero or not. The default value is \code{TRUE}.
#' @param zerolty single integer in [0,6] or one of the character strings \code{"blank",
#' "solid", "dashed", "dotted", "dotdash", "longdash", "twodash"}. Determines the line type of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{2}.
#' @param zerolwd single real, positive number. Determines the line width of the horizontal line at zero. Has no effect
#' if \code{drawzero = FALSE}. The default value is \code{1}.
#' @param zerocol single character string. Determines the color of the horizontal line at zero. Has no effect if \code{drawzero = FALSE}.
#' The default value is \code{"grey"}.
#' @param ... further arguments to be passed to \code{plot}.
#' @return Called for its side effects and returns invisibly.
#'
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' res <- shrinkTVP(y ~ x1 + x2, data)
#' plot(res$beta$beta_x1)
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.mcmc.tvp <- function(x, probs = c(0.025, 0.25, 0.75, 0.975),
                          shaded = TRUE, quantlines = FALSE,
                          shadecol = "skyblue", shadealpha = 0.5,
                          quantlty = 2, quantcol = "black", quantlwd = 0.5,
                          drawzero = TRUE, zerolty = 2, zerolwd = 1, zerocol = "grey", ...){

  # Input checking
  if (length(probs) == 0 || is.vector(probs) == FALSE | is.list(probs) == TRUE ){
    stop("probs has to be a vector")
  }

  # Check if all probs are specified correctly
  if (any(probs > 1) | any(probs < 0) | any(sapply(probs, numeric_input_bad_zer))){
    stop("all elements of probs have to be numbers between >= 0 and <= 1")
  }

  if (length(probs[probs < 0.5]) != length(probs[probs > 0.5])){
    stop ("There has to be an equal amount of probs above and below 0.5")
  }

  if (length(quantlines) == 0 || any(sapply(quantlines, bool_input_bad))){
    stop("all elements of quantlines have to be logical values")
  }

  if (length(shaded) == 0 || any(sapply(shaded, bool_input_bad))){
    stop("all elements of shaded have to be logical values")
  }

  if (is.character(shadecol)){
    if (any(sapply(shadecol, char_input_bad))){
      stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else if (is.numeric(shadecol)){
    if (any(sapply(shadecol, numeric_input_bad))){
      stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else {
    stop("shadecol has to be a string or a positive number or a vector of strings or positive numbers")
  }


  if (length(shadealpha) == 0 || any(shadealpha > 1) | any(shadealpha < 0) | any(sapply(shadealpha, numeric_input_bad_zer))){
    stop("all elements of shadealpha have to be numbers between >= 0 and <= 1")
  }

  if (length(quantlty) == 0 || any(sapply(quantlty, lty_input_bad))) {
    stop("every element of quantlty has to be an integer from 0 to 6 or one of 'blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', or 'twodash'")
  }

  if (is.character(quantcol)){
    if (any(sapply(quantcol, char_input_bad))){
      stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else if (is.numeric(quantcol)){
    if (any(sapply(quantcol, numeric_input_bad))){
      stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
    }
  } else {
    stop("quantcol has to be a string or a positive number or a vector of strings or positive numbers")
  }

  if (length(quantlwd) == 0 || any(sapply(quantlwd, numeric_input_bad))){
    stop("quantlwd has to be either a positive number or a string of positive numbers")
  }

  if (bool_input_bad(drawzero)){
    stop("drawzero has to be a logical")
  }

  if (lty_input_bad(zerolty)){
    stop("zerolty has to be an integer from 0 to 6 or one of 'blank', 'solid', 'dashed', 'dotted', 'dotdash', 'longdash', or 'twodash'")
  }

  if (numeric_input_bad(zerolwd)){
    stop("zerolwd has to be a positive number")
  }

  if (is.character(zerocol)){
    if (char_input_bad(zerocol)){
      stop("zerocol has to be a string or a positive number")
    }
  } else if (is.numeric(zerocol)){
    if (numeric_input_bad(zerocol)){
      stop("zerocol has to be a string or a positive number")
    }
  } else {
    stop("zerocol has to be a string or a positive number")
  }


  # Sort by ascending size
  probs <- c(probs, 0.5)
  probs <- sort(probs)

  # Calculate all quantiles
  if (length(probs) == 1){
    quants <- matrix()
  }

  # Create x_vec for plotting
  if (!is.null(attr(x, "index"))) {
    x_vec <- attr(x, "index")
  } else {
    x_vec <- 1:ncol(x)
  }


  startpoint <- ifelse(length(x_vec) == ncol(x), 1, 2)
  quants <- apply(x[, startpoint:ncol(x)], 2, quantile, probs)

  # Extract all user specified args and add standard ones that are missing
  args <- list(...)
  standard_args <- list(ylim = c(min(quants), max(quants)), type = "l", xlab = "", ylab = "")
  missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
  args[missing_args] <- standard_args[missing_args]

  # Add x and y to arguments
  args$x <- x_vec
  if (length(probs) == 1){
    args$y <- quants
  } else {
    args$y <- quants[median(1:nrow(quants)), ]
  }

  # Create plot
  do.call(plot, args)


  # Add horizontal line at 0
  if (drawzero){
    abline(h = 0, lty = zerolty, col = zerocol, lwd = zerolwd)
  }


  nint <- ((nrow(quants) - 1)/2)
  for (i in 1:nint){
    currcol <- ifelse(length(shadecol) > 1, rep_len(shadecol, nint)[i], shadecol)
    curralpha <- ifelse(length(shadealpha) > 1,  rep_len(shadealpha, nint)[i], shadealpha)

    currlty <-  ifelse(length(quantlty) > 1, rep_len(quantlty, nint)[i], quantlty)
    currbordcol <- ifelse(length(quantcol) > 1, rep_len(quantcol, nint)[i], quantcol)
    currlwd <-  ifelse(length(quantlwd) > 1, rep_len(quantlwd, nint)[i], quantlwd)
    currbord <- ifelse(length(quantlines) > 1, rep_len(quantlines, nint)[i], quantlines)

    polygon(c(x_vec, rev(x_vec)), c(quants[i, ], rev(quants[nrow(quants) + 1 - i, ])),
            col = adjustcolor(currcol, alpha.f = curralpha), border = FALSE)

    if (currbord == TRUE){
      lines(y = quants[i, ], x = x_vec, col = currbordcol, lwd = currlwd, lty = currlty)
      lines(y = quants[nrow(quants) + 1 - i, ], x = x_vec, col = currbordcol, lwd = currlwd, lty = currlty)
    }
  }

  # re paint line to be on top
  do.call(lines, args[names(args) %in% c("y", "x", "lwd", "lty", "col")])

}

#' Graphical summary of posterior distribution
#'
#' \code{plot.shrinkTVP} generates plots visualizing the posterior distribution.
#'
#' @param x a \code{shrinkTVP} object.
#' @param pars a character vector containing the names of the parameters to be visualized.
#' The names have to coincide with the names of the list elements of the \code{shrinkTVP}
#' object. Throws an error if any element of \code{pars} does not fulfill this criterium.
#' The default is \code{c("beta")}.
#' @param nplot positive integer that indicates the number of tvp plots to display on a single
#' page before a new page is generated. The default value is 3.
#' @param h_borders single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
#' the relative amount of space (the total amount summing up to 1) left blank on the left and right of the plot, in that order.
#' The default is \code{c(0.05, 0.05)}.
#' @param w_borders single real, positive number smaller than 0.5 or a vector containing two such numbers. Determines
#' the relative amount of space (the total amount summing up to 1) left blank at the top and bottom of the plot, in that order.
#' The default is \code{c(0.02, 0.02)}.
#' @param ... further arguments to be passed to the respective plotting functions.
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP(theta = c(0.2, 0, 0), beta_mean = c(1.5, -0.3, 0))
#' data <- sim$data
#'
#' output <- shrinkTVP(y ~ x1 + x2, data)
#' plot(output)
#' }
#'
#' ## Will produce an error because 'hello' is not a parameter in the model
#' \dontrun{
#' plot(output, pars = c("beta", "hello"))
#' }
#'
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkTVP <- function(x, pars = c("beta"), nplot = 3, h_borders = c(0.05, 0.05), w_borders = c(0.02, 0.02), ...){


  # Check that any pars are specified
  if (length(pars) == 0){
    stop("Please specify pars to plot.")
  }

  # Check if pars is of type "char"
  if (is.character(pars) == FALSE){
    stop("pars has to be of type characater")
  }

  # Check if all specified par values are allowed
  sample_names <- names(x)[sapply(x, function(y) attr(y, "type")) == "sample"]
  if (any(!pars %in% sample_names)){
    stop(paste0("pars can only contain the following values: ", paste(sample_names, collapse = ", ")))
  }

  # Check if nplot is okay
  if (int_input_bad(nplot)){
    stop("nplot has to be a single, positive integer")
  }

  if (nplot == 0){
    stop("nplot cant be 0")
  }

  if (any(sapply(h_borders, numeric_input_bad_zer)) | sum(h_borders) >= 1 | length(h_borders) > 2 |
      (length(h_borders) == 1 & h_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (any(sapply(w_borders, numeric_input_bad_zer)) | sum(w_borders) >= 1 | length(w_borders) > 2 |
      (length(w_borders) == 1 & w_borders > 0.5)[1]) {
    stop("h_borders has to be vector of length 2 containing numbers that sum up to less than 1 or a single number that is smaller than 0.5")
  }

  if (length(h_borders) == 1){
    h_borders <- rep(h_borders, 2)
  }

  if (length(w_borders) == 1){
    w_borders <- rep(w_borders, 2)
  }

  # Create sublist containing only parameters specified in par and copy attributes
  obj <- x[pars]
  mostattributes(obj) <- attributes(x)
  names(obj) <- names(x[pars])

  # Save previous user defined par attributes
  prev_par <- par(no.readonly = TRUE)

  # Idea: loop over all list elements and apply plot method specific to that object (mcmc or mcmc.tvp)
  for (i in 1:length(obj)){

    # TVP parameters will be a list, so we differentiate between TVP and non-TVP this way
    if (is.list(obj[[i]]) == TRUE){

      curr_nplot <- min(length(obj[[i]]), nplot)
      p_left <- length(obj[[i]])

      w_cut <- sum(w_borders)
      h_cut <- sum(h_borders)

      layout(rbind(0,
                   cbind(0, 1:curr_nplot, 0),
                   0),
             heights = c(h_borders[1], rep((1-h_cut)/curr_nplot, curr_nplot), h_borders[2]),
             widths = c(w_borders[1], 1 - w_cut, w_borders[2]))

      # Apply plotting method in a loop
      for (j in 1:length(obj[[i]])){

        if ((j-1) %% nplot == 0){
          # The first plot on a new page
          xaxt <- "n"
          curr_pos <- 0
        } else if ((j-1) %% nplot == (nplot - 1)){
          # The last plot on a page
          xaxt <- "s"
        } else if (p_left < nplot & j == length(obj[[i]])){
          # The last plot on the last page
          xaxt <- "s"
        } else {
          # Any plot in between two others
          xaxt <- "n"
        }
        par(mar = c(0, 2.6, 0, 0), mgp = c(1.6, .6, 0))

        # Extract all user specified args
        args <- list(...)
        # Augment user specified arguments with standard ones
        standard_args <- list(ylab = paste(names(obj)[i], "of", attr(obj, "colnames")[j]), xaxt = xaxt)
        missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
        args[missing_args] <- standard_args[missing_args]

        # Add object to be plotted to args list
        args$x <- obj[[i]][[j]]

        # Plot
        do.call(plot, args)

        # Wait for user to move on
        if (j %% nplot == 0 & j != length(obj[[i]])){
          readline("Hit <Return> to see next plot: ")
        }

        p_left <- p_left - j * nplot

      }

      # Move to next series of plots if current parameter has been thoroughly plotted
      if (i < length(obj)){
        readline("Hit <Return> to see next plot: ")
        par(prev_par)
      }

    } else {

      # In non-TVP case, simply use coda's plot method
      plot(obj[[i]], ...)

      if (i < length(obj)){
        readline("Hit <Return> to see next plot: ")
      }
    }
  }

  # Restore user defined values of par
  par(prev_par)
}

#' Graphical summary of posterior predictive density
#'
#' \code{plot.shrinkTVP_forc} generates plots visualizing the posterior predictive density generated by \code{forecast_shrinkTVP}.
#'
#' @param x a \code{shrinkTVP_forc} object.
#' @param showgap if \code{showgap = FALSE}, the gap between the historical observations and the forecasts is removed.
#' The default value is \code{FALSE}.
#' @param ... further arguments to be passed to \code{plot}.
#'
#' @return Called for its side effects and returns invisibly.
#' @examples
#' \donttest{
#' set.seed(123)
#' sim <- simTVP()
#'
#' train <- sim$data[1:190, ]
#' test <- sim$data[191:200, ]
#'
#' res <- shrinkTVP(y ~ x1 + x2, train)
#'
#' forecast <- forecast_shrinkTVP(res, test)
#' plot(forecast)
#' lines(sim$data$y, col = "forestgreen")
#'
#' }
#' @author Peter Knaus \email{peter.knaus@@wu.ac.at}
#' @family plotting functions
#' @export
plot.shrinkTVP_forc<- function(x, showgap = FALSE, ...){

  args <- list(...)

  if ("probs" %in% names(args)){
    probs <- args$probs

    # Input checking
    if (length(probs) == 0 || is.vector(probs) == FALSE | is.list(probs) == TRUE ){
      stop("probs has to be a vector")
    }

    # Check if all probs are specified correctly
    if (any(probs > 1) | any(probs < 0) | any(sapply(probs, numeric_input_bad_zer))){
      stop("all elements of probs have to be numbers between >= 0 and <= 1")
    }

    if (length(probs[probs < 0.5]) != length(probs[probs > 0.5])){
      stop ("There has to be an equal amount of probs above and below 0.5")
    }
  } else {
    probs <- c(0.025, 0.25, 0.75, 0.975)
  }

  prelim_probs <- range(probs)
  prelim_quants <- apply(x$y_pred, 2, quantile, prelim_probs)

  index_old <- attr(x, "index")
  spacing <- diff(index_old)[1]

  index_new <- seq(from = index_old[length(index_old)] + spacing,
                   to = index_old[length(index_old)] + ncol(x$y_pred) * spacing,
                   by = spacing)

  stand_ylim <- c(min(prelim_quants, x$y_orig), max(prelim_quants, x$y_orig))
  standard_args <- list(ylab = "", xlim = c(index_old[1], index_new[length(index_new)]),
                        ylim = stand_ylim)
  missing_args <- names(standard_args)[!names(standard_args) %in% names(args)]
  args[missing_args] <- standard_args[missing_args]

  pred <- x$y_pred
  if (showgap == FALSE){
    index_new <- c(index_new[1] - spacing, index_new)
    pred <- cbind(x$y_orig[length(x$y_orig)], pred, row.names = NULL)
  }

  args$x <- pred
  attr(args$x, "index") <- index_new

  # Plot
  do.call(plot.mcmc.tvp, args)
  args$y <- x$y_orig
  args$x <- index_old
  args$probs <- NULL
  do.call(lines, args)

}
