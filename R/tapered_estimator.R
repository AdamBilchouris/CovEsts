#' Compute Normalisation Factor
#'
#' This helper function is used in the computation of the normalisation factor the function [tapered_est],
#' \deqn{H_{2, n}(0) = \sum_{s=1}^{n} a((s - 1/2) / n; \rho)^{2}, }
#' where \eqn{a(\cdot; \cdot)} is a window function.
#'
#' @param n The sample size.
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window_ec] function to be used. Possible values are:
#' tukey, triangular, power_sine, blackman_window, hann_poisson, welch. Alternatively, a custom window function can be provided, see the example for [taper].
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#'
#' @return A single value being \eqn{H_{2, n}(0)}.
#'
#' @examples
#' \dontrun{
#' H2n(3, 0.6, "tukey")
#' }
H2n <- function(n, rho, window_name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(n), n >= 1, n %% 1 == 0, is.numeric(rho), rho > 0, rho <= 1, length(rho) == 1, is.logical(custom_window))
  window_name <- match.arg(window_name)
  sSeq <- 1:n
  hSeq <- taper(((sSeq - 1/2) / n), rho, window_name, window_params, custom_window)
  return(sum(hSeq))
}

#' Compute the Function \eqn{a(x; \rho).}
#'
#' This helper function computes the taper function for a given window function as
#' \deqn{a(x; \rho) = \left\{
#' \begin{array}{ll}
#'    w(2x/\rho)     & 0 \leq x < \frac{1}{2} \rho, \\
#'    1              & \frac{1}{2}\rho \leq x \leq \frac{1}{2} \\
#'    a(1 - x; \rho) & \frac{1}{2} < x \leq 1
#' \end{array} ,
#' \right. }
#' where \eqn{w(\cdot)} is a continuous increasing function with \eqn{w(0)=0, w(1)=1,}
#' \eqn{\rho \in (0, 1],} and \eqn{x \in [0, 1].} The possible window function choices are found in [window_ec].
#'
#' @param x A vector of numbers between 0 and 1 (inclusive).
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window_ec] function to be used. Possible values are:
#' tukey, triangular, power_sine, blackman_window, hann_poisson, welch. Alternatively, a custom window function can be provided, see the example.
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#'
#' @return A vector of taper values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.2, 0.3)
#' taper(X, 0.5, "tukey")
#' curve(taper(x, 1, "tukey"), from = 0, to = 1)
#' curve(taper(x, 1, "power_sine", c(4)), from = 0, to = 1)
#' my_taper <- function(x, ...) {
#'   return(x)
#' }
#' taper(X, 0.5, my_taper, custom_window = TRUE)
taper <- function(x, rho, window_name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), is.numeric(rho),
            length(rho) == 1, rho > 0, rho <= 1, is.logical(custom_window))
  # after 1/2 , we have a(1-x), treat it as that.
  x_new <- pmin(x, 1 - x)

  # Check which are not in [rho/2, 1/2] : values in this are mapped to 1.
  x_new_idx <- x_new < (rho / 2)

  # allocate 1 as we return it for non in x_new_idx
  retVec <- rep(1, length(x))

  if(custom_window) {
    stopifnot(is.function(window_name))
    retVec[x_new_idx] <- window_name((2 * x_new[x_new_idx]) / rho, window_params)
  }

  else {
    window_name <- match.arg(window_name)
    retVec[x_new_idx] <- window_ec((2 * x_new[x_new_idx]) / rho, window_name, c(window_params[1]))
  }
  return(retVec)
}

#' Compute the Estimated Tapered Autocovariance Function over a Set of Lags.
#'
#' This function computes the tapered autocovariance over a set of lags.
#'
#' @details
#' This function computes the estimated tapered autocovariance over a set of lags,
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n}(0))^{-1} \sum_{j=1}^{N-h}  (X(j) - \bar{X} ) ( X(j + h) - \bar{X} ) a((j - 1/2) / N; \rho) a((j + h - 1/2) / N; \rho) ,}
#' where \eqn{a(\cdot)} is a window function, \eqn{\rho \in (0, 1]}  is a scale parameter.
#' This estimator takes into account the edge effect during estimation, assigning a lower weight to values closer to the boundaries and higher weights for observations closer to the middle.
#' This estimator is positive-definite and asymptotically unbiased.
#'
#' @references
#' Dahlhaus R. & Künsch, H. (1987). Edge Effects and Efficient Parameter Estimation for Stationary Random Fields. Biometrika 74(4), 877-882. 10.1093/biomet/74.4.877
#'
#' @param X A vector representing observed values of the time series.
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window_ec] function to be used. Possible values are:
#' tukey, triangular, sine, power_sine, blackman_window, hann_poisson, welch. Alternatively, a custom window_ec function can be provided, see the example in [taper].
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param x A vector of lag indices. Defaults to the sequence \code{0:length(X)}. Must be at least as large as \code{maxLag + 1}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#'
#' @return A \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the autocovariance/autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the \code{type} parameter.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'tapered_est'.}
#' }
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' tapered_est(X, 0.5, "tukey", maxLag = 2)
tapered_est <- function(X, rho, window_name = c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"), window_params = c(1), maxLag = length(X) - 1, x = 0:length(X), type = c("autocovariance", "autocorrelation"), meanX = mean(X), custom_window = FALSE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(rho), length(rho) == 1, rho > 0, rho <= 1,
            is.logical(custom_window),  length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            is.numeric(x), all(!is.na(x)), length(x) > 0)

  type <- match.arg(type)

  h2n <- H2n(length(X), rho, window_name, window_params, custom_window)
  covVals <- rep(0, maxLag + 1)

  n <- length(X)
  taperVals_t <- taper(((1:n) - 1/2)/n, rho, window_name, window_params, custom_window)

  # auxiliary variable, compute the ACF of this new sequence instead.
  Y <- (X - mean(X)) * taperVals_t

  # compute the PD estimated ACF of Y, rescale to 1/h2n
  covVals <- standard_est(Y, maxLag = maxLag, type = type, meanX = 0)$acf * n / h2n

  if(type == 'autocorrelation') {
    covVals <- covVals / covVals[1]
  }

  res <- list(acf = covVals, lags = x[1:(maxLag + 1)], est_type = type, est_used = 'tapered_est')
  return(structure(res, class = "CovEsts"))
}
