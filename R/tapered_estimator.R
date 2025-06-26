#' Compute Taper
#'
#' This helper function computes the taper function for a given window.
#' \deqn{a(x; \rho) = \left\{
#' \begin{array}{ll}
#'    w(2x/\rho)     & 0 \leq x < \frac{1}{2} \rho, \\
#'    1              & \frac{1}{2}\rho \leq x \leq \frac{1}{2} \\
#'    a(1 - x; \rho) & \frac{1}{2} < x \leq 1
#' \end{array} ,
#' \right. }
#' where \eqn{w(\cdot)} is a continuous increasing function with \eqn{w(0)=0, w(1)=1,}
#' \eqn{\rho \in (0, 1],} and \eqn{x \in [0, 1].} The possible window choices are found in [window].
#'
#' @param x A number between 0 and 1 (inclusive).
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window] function to be used. Possible values are:
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch". Alternatively, a custom window function can be provided, see the example.
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#'
#' @return A single value being the taper function value at x.
#' @export
#'
#' @examples
#' x <- 0.4
#' taper_single(x, 0.5, "tukey")
#' my_taper <- function(x, ...) {
#'   return(x)
#' }
#' taper_single(x, 0.5, "my_taper", custom_window = TRUE)
taper_single <- function(x, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) == 1, x >= 0, x <= 1, is.numeric(rho), rho > 0, rho <= 1, length(rho) == 1, is.logical(custom_window))

  if(custom_window) {
    if(x >= 0 && x < ((1/2) * rho)) {
      return(get(window_name)((2*x) / rho, window_params))
    }
    else if(x >= ((1/2) * rho) && x <= (1/2)) {
      return(1)
    }
    else if(x > (1/2) && x <= 1) {
      x <- 1 - x
      if(x >= 0 && x < ((1/2) * rho)) {
        return(get(window_name)((2*x) / rho, window_params))
      }
      else if(x >= ((1/2) * rho) && x <= (1/2)) {
        return(1)
      }
    }
  }

  else {
    stopifnot(window_name %in% c("tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch"))
    if(x >= 0 && x < ((1/2) * rho)) {
      return(get("window")((2*x) / rho, window_name, c(window_params[1])))
    }
    else if(x >= ((1/2) * rho) && x <= (1/2)) {
      return(1)
    }
    else if(x > (1/2) && x <= 1) {
      x <- 1 - x
      if(x >= 0 && x < ((1/2) * rho)) {
        return(get("window")((2*x) / rho, window_name, c(window_params[1])))
      }
      else if(x >= ((1/2) * rho) && x <= (1/2)) {
        return(1)
      }
    }
  }
  return(NA)
}

#' Compute \eqn{H{2, n}}.
#'
#' This helper function is used in the computation of the function [tapered_cov_single].
#' \deqn{H_{2, n}(0) = \sum_{s=1}^{n} a((s - 1/2) / n; \rho)^{2}, }
#' where \eqn{a(\cdot; \cdot)} is a window function.
#'
#' @param n The sample size.
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window] function to be used. Possible values are:
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch". Alternatively, a custom window function can be provided, see the example.
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#'
#' @return A single value being \eqn{H_{2, n}(0)}.
#' @export
#'
#' @examples
#' H2n(3, 0.6, "tukey")
H2n <- function(n, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(n), n >= 1, n %% 1 == 0, is.numeric(rho), rho > 0, rho <= 1, length(rho) == 1, is.logical(custom_window))
  sSeq <- 1:n
  hSeq <- taper(((sSeq - 1/2) / n), rho, window_name, window_params, custom_window)
  return(sum(hSeq))
}

#' Compute the Function \eqn{a(x; \rho).}
#'
#' This function repeatedly calls [taper_single] (refer to its manual).
#'
#' @param x A vector of number between 0 and 1 (inclusive).
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window] function to be used. Possible values are:
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch". Alternatively, a custom window function can be provided, see the example.
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
taper <- function(x, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), is.numeric(rho),
            length(rho) == 1, rho > 0, rho <= 1, is.logical(custom_window))
  retTaper <- c()
  for(i in 1:length(x)) {
    retTaper <- c(retTaper, taper_single(x[i], rho, window_name, window_params, custom_window))
  }
  return(retTaper)
}

#' Computes the Tapered Autocovariance for a Single Lag.
#'
#' This helper function computes the tapered autocovariance for a single lag \eqn{h},
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n})^{-1} \sum_{t=1}^{N-h}  (X(t_{i}) - \bar{X} ) ( X(t_{i} + h) - \bar{X} ) [ a((t_{i} - 1/2) / n; \rho) a((t_{i} + h - 1/2) / n; \rho)  ] ,}
#' where \eqn{a(\cdot)} is a window function, \eqn{\rho} is a scale parameter. The taper functions are precomputed in [compute_tapered_est].
#'
#' @param X A vector representing observed values of the time series.
#' @param meanX The average value of the \code{X}.
#' @param h The lag at which the tapered autocovariance function is computed at.
#' @param h2n The value of \eqn{H_{2, n}(0)}, computed within [compute_tapered_est].
#' @param taperVals_t The taper values for each index of the process, computed within [compute_tapered_est].
#' @param taperVals_h The taper values shifted by the lag, computed within [compute_tapered_est].
#'
#' @return The tapered autocovariance function at the specified lag.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' tapered_cov_single(X, mean(X), 9, 2.5, c(0.75, 1, 0.75), c(0.75, 1, 0.75))
tapered_cov_single <- function(X, meanX, h, h2n, taperVals_t, taperVals_h) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(meanX), length(meanX) == 1,
            is.numeric(h), length(h) == 1, h %% 1 == 0, is.numeric(h2n), length(h2n) == 1,
            is.numeric(taperVals_t), length(taperVals_t) >= 1, all(abs(taperVals_t) <= 1), all(abs(taperVals_t) >= 0),
            is.numeric(taperVals_h), length(taperVals_h) >= 1, all(abs(taperVals_h) <= 1), all(abs(taperVals_h) >= 0))

  covVals <- c()
  n <- length(X)

  for(i in 1:(n - h)) {
    covVals <- c(covVals, (X[i] - meanX) * (X[i + h] - meanX) * taperVals_t[i] * taperVals_h[i])
  }
  return(sum(covVals) / h2n)
}

#' Compute the Estimated Tapered Autocovariance Function over a Set of Lags.
#'
#' This function computes the tapered autocovariance over a set of lags,
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n}(0))^{-1} \sum_{j=1}^{N-h}  (X(j) - \bar{X} ) ( X(j + h) - \bar{X} ) a((j - 1/2) / N; \rho) a((j + h - 1/2) / N; \rho) ,}
#' where \eqn{a(\cdot)} is a window function, \eqn{\rho \in (0, 1]} is a scale parameter.
#' For each lag, the tapered autocovariance is computed using the function [tapered_cov_single].
#'
#' @details
#' This function computes the estimated tapered autocovariance over a set of lags,
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n}(0))^{-1} \sum_{j=1}^{N-h}  (X(j) - \bar{X} ) ( X(j + h) - \bar{X} ) a((j - 1/2) / N; \rho) a((j + h - 1/2) / N; \rho) ,}
#' where \eqn{a(\cdot)} is a window function, \eqn{\rho \in (0, 1]}  is a scale parameter.
#' This estimator considers the edge effect during estimation, assigning a lower weight to values closer to the boundaries and higher weights for observations closer to the middle.
#' This estimator is positive-definite and asymptotically unbiased.
#'
#' Internally, this function calls [tapered_cov_single] for each lag \eqn{h}.
#'
#' @references
#' Dahlhaus R., Künsch, H. (1987). Edge Effects and Efficient Parameter Estimation for Stationary Random Fields. In Biometrika (Vol. 74 Issue 4, pp. 877-882). 10.1093/biomet/74.4.877
#'
#' @param X A vector representing observed values of the time series.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param rho A scale parameter in \eqn{(0, 1].}
#' @param window_name The name of the [window] function to be used. Possible values are:
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch". Alternatively, a custom window function can be provided, see the example.
#' @param window_params A vector of parameters of the window function.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#'
#' @return  A vector whose values are the estimated tapered autocovariances.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_tapered_est(X, 2, 0.5, "tukey")
compute_tapered_est <- function(X, maxLag, rho, window_name, window_params = c(1), custom_window = FALSE, type = 'autocovariance',  meanX = mean(X)) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(rho), length(rho) == 1, rho > 0, rho <= 1,
            is.logical(custom_window),  length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            type %in% c('autocovariance', 'autocorrelation'))

  h2n <- H2n(length(X), rho, window_name, window_params, custom_window)
  covVals <- c()

  n <- length(X)
  taperVals_t <- taper(((1:n) - 1/2)/n, rho, window_name, window_params, custom_window)

  for(i in 0:maxLag) {
    taperVals_h <- taperVals_t[(1 + i):length(taperVals_t)]

    covVals <- c(covVals, tapered_cov_single(X, meanX, i, h2n, taperVals_t, taperVals_h))
  }

  if(type == 'autocorrelation') {
    return(covVals / covVals[1])
  }

  return(covVals)
}
