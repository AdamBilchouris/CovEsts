#' Plot Method for CovEsts Objects
#'
#' @description
#' This function plots a CovEsts object, over the lags used to estimate the autocovariance values.
#'
#' @param x A CovEst S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param xlab Defaults to \code{'Lag (h)'}, see [graphics::title] for more information.
#' @param ylab Defaults to \code{'ACF'}, see [graphics::title] for more information.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A plot of a CovEsts S3 object.
#' @export
#'
#' @examples
#' plot(standard_est(c(1, 2, 3)))
plot.CovEsts <- function(x, type = 'l', xlab = 'Lag (h)', ylab = 'ACF', ...) {
  # add logic to see if it is autocorrelation, autocovariance, partial
  # change plot y-axis accordingly.
  plot(x = x$lags, y = x$acf, type = type, xlab = xlab, ylab = ylab, ...)
}

#' Lines Method for CovEsts Objects
#'
#' @description
#' This function plots a CovEsts object as a line onto another plot.
#'
#' @param x A CovEst S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A line of a CovEsts S3 object.
#' @export
#'
#' @examples
#' plot(standard_est(c(1, 2, 3)))
#' lines(standard_est(c(1, 2, 3), pd = FALSE))
lines.CovEsts <- function(x, type = "l", ...) {
  lines(x = x$lags, y = x$acf,  type = type, ...)
}

#' Print Method for CovEsts Objects
#'
#' @description This function plots a CovEsts object, printing the estimator
#' type, estimator used, lags and values.
#'
#' @param x A CovEst S3 object.
#' @param n_head The number of 'head' values to print, defaults to 5.
#' @param n_tail The number of 'tail' values to print, defaults to 5.
#' @param digits The number of digits to print. Defaults to \code{NULL}. Note, this value is passed into [base::format].
#' @param ... Addition printing parameters, see [base::print].
#'
#' @return Prints a CovEsts object.
#' @export
#'
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @examples
#' print(standard_est(c(1, 2, 3)))
print.CovEsts <- function(x, n_head = 5, n_tail = 5, digits = NULL, ...) {
  cat("\n Estimated (Partial) Autocorrelation/Autocovariance\n\n")
  cat("Type:      ", x$est_type, "\n")
  cat("Estimator: ", x$est_used, "\n")

  full_df <- format(data.frame(Lag = x$lags, Estimate = x$acf), digits = digits)

  if(length(x$acf) <= (n_head + n_tail)) {
    print(full_df, row.names = FALSE, right = FALSE)
  }

  else {
    head_df <- utils::head(full_df, n_head)
    dots_df <- data.frame(Lag = "...", Estimate = "...")
    tail_df <- utils::tail(full_df, n_tail)

    print(rbind(head_df, dots_df, tail_df), row.names = FALSE, right = FALSE)
  }

  cat("\n")
  invisible(x)
}

#' Plot Method for VarioEsts Objects
#'
#' @description
#' This function plots a VarioEsts object, over the lags used to estimate the variogram values.
#'
#' @param x A VarioEsts S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param xlab Defaults to \code{'Lag (h)'}, see [graphics::title] for more information.
#' @param ylab Defaults to \code{expression(gamma(h))}, see [graphics::title] for more information.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A plot of a VarioEsts S3 object.
#' @export
#'
#' @examples
#' plot(to_vario(standard_est(c(1, 2, 3))))
plot.VarioEsts <- function(x, type = 'l', xlab = 'Lag (h)', ylab = expression(gamma(h)), ...) {
  plot(x = x$lags, y = x$vario, type = type, xlab = xlab, ylab = ylab, ...)
}

#' Lines Method for VarioEsts Objects
#'
#' @description
#' This function plots a VarioEsts object as a line onto another plot.
#'
#' @param x A VarioEsts S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A line of a VarioEsts S3 object.
#' @export
#'
#' @examples
#' plot(to_vario(standard_est(c(1, 2, 3))))
#' lines(to_vario(standard_est(c(1, 2, 3), pd = FALSE)))
lines.VarioEsts <- function(x, type = "l", ...) {
  lines(x = x$lags, y = x$vario, type = type, ...)
}

#' Print Method for VarioEsts Objects
#'
#' @description This function plots a VarioEsts object, printing the lags and values.
#'
#' @param x A VarioEsts S3 object.
#' @param n_head The number of 'head' values to print, defaults to 5.
#' @param n_tail The number of 'tail' values to print, defaults to 5.
#' @param digits The number of digits to print. Defaults to \code{NULL}. Note, this value is passed into [base::format].
#' @param ... Addition printing parameters, see [base::print].
#'
#' @return Prints a VarioEsts object.
#' @export
#'
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @examples
#' print(to_vario(standard_est(c(1, 2, 3))))
print.VarioEsts <- function(x, n_head = 5, n_tail = 5, digits = NULL, ...) {
  cat("\n Estimated Variogram\n\n")

  full_df <- format(data.frame(Lag = x$lags, Estimate = x$vario), digits = digits)

  if(length(x$vario) <= (n_head + n_tail)) {
    print(full_df, row.names = FALSE, right = FALSE)
  }

  else {
    head_df <- utils::head(full_df, n_head)
    dots_df <- data.frame(Lag = "...", Estimate = "...")
    tail_df <- utils::tail(full_df, n_tail)

    print(rbind(head_df, dots_df, tail_df), row.names = FALSE, right = FALSE)
  }

  cat("\n")
  invisible(x)
}

#' Plot Method for BootEsts Objects
#'
#' @description
#' This function plots a BootEsts object. It plots the confidence region, average bootstrap estimate, and the original nonbootstrapped estimate.
#'
#' @param x A BootEsts S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param xlab Defaults to \code{'Lag (h)'}, see [graphics::title] for more information.
#' @param ylab Defaults to \code{'ACF'}, see [graphics::title] for more information.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A plot of a BootEsts S3 object.
#' @export
#' @importFrom grDevices rgb
#'
#' @examples
#' plot(block_bootstrap(c(1, 2, 3), 2))
plot.BootEsts <- function(x, type = 'l', xlab = 'Lag (h)', ylab = 'ACF', ...) {
  plot(x$lags, x$acf_orig, type = "n", xlab = 'Lag (h)', ylab = 'ACF', ...)
  polygon(c(rev(x$lags), x$lags), c(rev(x$conf_upper), x$conf_lower), col = rgb(0.9, 0.9, 0.9, 0.6))
  lines(x$lags, x$acf_orig, type = "l", lwd = 2)
  lines(x$lags, x$acf_avg, lty = 2, col = 2, lwd = 2)

  graphics::legend("topright", inset=0, title=paste0(toupper(substring(x$boot_type, 1, 1)), substring(x$boot_type, 2), " Block Bootstrap"),
                   c("Original Estimate", "Bootstrapped Average", paste0((1 - x$alpha) * 100, "% Confidence Region")),
                   col=c(1, 2, rgb(0.9, 0.9, 0.9, 0.6)), lwd=c(2, 2, NA), lty=c(1, 2, NA), density=c(0, 0, NA), fill=c(NA, NA, "#AAAAAAAA"), border=c(NA, NA, "#000000FF"),
                   y.intersp=1, cex=1)
}

#' Lines Method for BootEsts Objects
#'
#' @description
#' This function plots a BootEsts object as a line onto another plot.
#'
#' @param x A BootEsts S3 object.
#' @param type Defaults to \code{'l'}, see [base::plot] for the possible values.
#' @param ... Additional plotting arguments, refer to [graphics::par].
#'
#' @return A line of a BootEsts S3 object.
#' @export
#'
#' @examples
#' plot(block_bootstrap(c(1, 2, 3), 2))
#' lines(block_bootstrap(c(1, 2, 3), 2, pd = FALSE))
lines.BootEsts  <- function(x, type = "l", ...) {
  lines(x = x$lags, y = x$acf_avg, type = type, ...)
}

#' Print Method for BootEsts Objects
#'
#' @description This function plots a BootEsts object, printing estimtor type, bootstrap type, estimator used, the bootstrap average estimate, the original estimate,
#' and the pointwise confidence interval.
#'
#' @param x A BootEsts S3 object.
#' @param n_head The number of 'head' values to print, defaults to 5.
#' @param n_tail The number of 'tail' values to print, defaults to 5.
#' @param digits The number of digits to print. Defaults to \code{NULL}. Note, this value is passed into [base::format].
#' @param ... Addition printing parameters, see [base::print].
#'
#' @return Prints a BootEsts object.
#' @export
#'
#' @importFrom utils head
#' @importFrom utils tail
#'
#' @examples
#' print(block_bootstrap(c(1, 2, 3), 2))
print.BootEsts <- function(x, n_head = 5, n_tail = 5, digits = NULL, ...) {
  cat("\n Estimated Bootstrap Autocorrelation/Autocovariance\n\n")
  cat("Type:      ", x$est_type, "\n")
  cat("Boot Type: ", x$boot_type, "\n")
  cat("Estimator: ", x$est_used, "\n")

  avg_df <- format(data.frame(Lag = x$lags, Estimate = x$acf_avg), digits = digits)
  orig_df <- format(data.frame(Lag = x$lags, Estimate = x$acf_orig), digits = digits)
  conf_df <- format(data.frame(Lag = x$lags, Conf_Lower = x$conf_lower, Conf_Upper = x$conf_upper), digits = digits)

  if(length(x$acf_avg) <= (n_head + n_tail)) {
    cat("Average ACF boostrap estimate\n")
    print(avg_df, row.names = FALSE, right = FALSE)
    cat("\nOriginal ACF estimate\n")
    print(orig_df, row.names = FALSE, right = FALSE)
    cat("\nBootstrap pointwise 1 - ", x$alpha, "% confidence interval\n")
    print(conf_df, row.names = FALSE, right = FALSE)
  }

  else {
    cat("Average ACF boostrap estimate\n")
    head_df <- utils::head(avg_df, n_head)
    dots_df <- data.frame(Lag = "...", Estimate = "...")
    tail_df <- utils::tail(avg_df, n_tail)
    print(rbind(head_df, dots_df, tail_df), row.names = FALSE, right = FALSE)

    cat("\nOriginal ACF estimate\n")
    head_df <- utils::head(orig_df, n_head)
    dots_df <- data.frame(Lag = "...", Estimate = "...")
    tail_df <- utils::tail(orig_df, n_tail)
    print(rbind(head_df, dots_df, tail_df), row.names = FALSE, right = FALSE)

    cat("\nBootstrap pointwise 1 - ", x$alpha, "% confidence interval\n", sep = "")
    head_df <- utils::head(conf_df, n_head)
    dots_df <- data.frame(Lag = "...", Conf_Lower = "...", Conf_Upper = "...")
    tail_df <- utils::tail(conf_df, n_tail)
    print(rbind(head_df, dots_df, tail_df), row.names = FALSE, right = FALSE)
  }
  cat("\n")
  invisible(x)
}


#' as.double Method for CovEsts Objects
#'
#' @description
#' This function provides a method for as.double for CovEsts objects, allowing for as.numeric to be called.
#' It extracts the acf values from the object.
#'
#' @param x A CovEsts S3 object.
#' @param ... Refer to [base::double].
#'
#' @return The acf values from the CovEsts object.
#' @export
#'
#' @examples
#' as.numeric(standard_est(c(1, 2, 3)))
as.double.CovEsts <- function(x, ...) {
  return(x$acf)
}
