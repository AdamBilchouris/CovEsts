#' Compute Taper
#'
#' This computes the taper function for a given window.
#' \deqn{a(x; \rho) =
#' \begin{cases}
#'    w(xu/\rho)     & 0 \leq x < \frac{1}{2} \rho, \\
#'    1              & \frac{1}{2}\rho \leq x \leq \frac{1}{2} \\
#'    a(1 - x; \rho) & \frac{1}{2} < x \leq 1
#' \end{cases} ,
#' }
#' where \eqn{w(\cdot)} is a continuous increasing function with \eqn{w(0)=0, w(1)=1,}
#' \eqn{\rho > 0,} and \eqn{x \in [0, 1].}
#'
#' @param x A number between 0 and 1 (inclusive).
#' @param rho A scale paramter, greater than 0.
#' @param window_name The window function to use, one of
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch".
#' @param window_params Any parameters of the window function.
#' @param custom_window Whether or not a custom window is used (boolean). Put its parameters in window_params.
#'
#' @return A single value being the taper function value at x.
#' @export
#'
#' @examples
#' x <- 0.4
#' compute_taper_single(x, 0.5, "tukey")
#' my_taper <- function(x, ...) {
#'   return(x)
#' }
#' compute_taper_single(x, 0.5, "my_taper", custom_window = TRUE)
compute_taper_single <- function(x, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) == 1, x >= 0, x <= 1, is.numeric(rho), rho > 0)

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

#' Compute \eqn{H_{2, n}}.
#'
#' This is used in the computation of [compute_tapered_cov_single].
#' \deqn{H_{2, n}(0) = \sum_{s=1}^{n} a((s - 1/2) / n)^2 , }
#' where \eqn{a} is a window function.
#'
#' @param n The number of samples.
#' @param rho A scale paramter, greater than 0.
#' @param window_name The window function to use, one of
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch".
#' @param window_params Any parameters of the window function.
#' @param custom_window Whether or not a custom window is used (boolean). Put its parameters in window_params.
#'
#' @return A single value being \eqn{H_{2, n}(0)}.
#' @export
#'
#' @examples
#' compute_h2n(3, 0.6, "tukey")
compute_h2n <- function(n, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(n), n >= 1)
  sSeq <- 1:n
  hSeq <- compute_taper(((sSeq - 1/2) / n), rho, window_name, window_params, custom_window)
  return(sum(hSeq))
}

#' Compute the function \eqn{a(x; \rho).}
#'
#' This function repeatedly calls [compute_taper_single].
#'
#' @param x A vector of number between 0 and 1 (inclusive).
#' @param rho A scale paramter, greater than 0.
#' @param window_name The window function to use, one of
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch".
#' @param window_params Any parameters of the window function.
#' @param custom_window Whether or not a custom window is used (boolean). Put its parameters in window_params.
#'
#' @return A vector of taper values.
#' @export
#'
#' @examples
#' X <- c(0.1, 0.2, 0.3)
#' compute_taper(X, 0.5, "tukey")
compute_taper <- function(x, rho, window_name, window_params=c(1), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), is.numeric(rho),
            length(rho) == 1, is.logical(custom_window))
  retTaper <- c()
  for(i in 1:length(x)) {
    retTaper <- c(retTaper, compute_taper_single(x[i], rho, window_name, window_params, custom_window))
  }
  return(retTaper)
}

#' Computes the tapered covariance for a single lag.
#'
#' This computes the tapered covariance for a single \eqn{h},
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n})^{-1} \sum_{t=1}^{N-h}  (X(t_{i}) - \bar{X} ) ( X(t_{i} + h) - \bar{X} ) [ a((t_{i} - 1/2) / n; \rho) a((t_{i} + h - 1/2) / n; \rho)  ] ,}
#' where \eqn{a} is a window function, \eqn{\rho} is a scale parameter. The taper functions are precomputed in [compute_tapered_cov].
#'
#' @param X A numeric vector representing the process.
#' @param meanX The mean of the process.
#' @param h The lag at which the tapered covariance function is computed at.
#' @param h2n The value of \eqn{H_{2, n}}, computed within [compute_tapered_cov].
#' @param taperVals_t The taper values for each index of the process, computed within [compute_tapered_cov].
#' @param taperVals_h The taper values shifted by the lag, computed within [compute_tapered_cov].
#'
#' @return The tapered covariance function at that lag.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_tapered_cov_single(X, mean(X), 9, 2.5, c(0.75, 1, 0.75), c(0.75, 1, 0.75))
compute_tapered_cov_single <- function(X, meanX, h, h2n, taperVals_t, taperVals_h) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(meanX), length(meanX) == 1,
            is.numeric(h), length(h) == 1, is.numeric(h2n), length(h2n) == 1,
            is.numeric(taperVals_t), length(taperVals_t) >= 1, all(abs(taperVals_t) <= 1), all(abs(taperVals_t) >= 0),
            is.numeric(taperVals_h), length(taperVals_h) >= 1, all(abs(taperVals_h) <= 1), all(abs(taperVals_h) >= 0))

  covVals <- c()
  n <- length(X)

  for(i in 1:(n - h)) {
    covVals <- c(covVals, (X[i] - meanX) * (X[i + h] - meanX) * taperVals_t[i] * taperVals_h[i])
  }
  return(sum(covVals) / h2n)
}

#' Compute the tapered covariance function over a set of lags.
#'
#' This computes the tapered covariance over a set of lags,
#' \deqn{\widehat{C}_{N}^{a} (h) = (H_{2, n})^{-1} \sum_{t=1}^{N-h}  (X(t_{i}) - \bar{X} ) ( X(t_{i} + h) - \bar{X} ) [ a((t_{i} - 1/2) / n; \rho) a((t_{i} + h - 1/2) / n; \rho)  ] ,}
#' where \eqn{a} is a window function, \eqn{\rho} is a scale parameter.
#' For each lag, the tapered covariance is computed in [compute_tapered_cov_single].
#'
#' @param X A numeric vector representing the process.
#' @param maxLag The maximum lag at which the covariance function is to be computed. Must be less than the length of the process.
#' @param rho A scale parameter greater than zero.
#' @param window_name The window function to use, one of
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch".
#' @param window_params Any parameters of the window function.
#' @param custom_window Whether or not a custom window is used (boolean). Put its parameters in window_param
#'
#' @return A numeric vector representing values of the tapered covariance function.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_tapered_cov(X, 2, 0.5, "tukey")
compute_tapered_cov <- function(X, maxLag, rho, window_name, window_params = c(1), custom_window = FALSE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1, maxLag <= (length(X) - 1),
            is.numeric(rho), length(rho) == 1, is.logical(custom_window))

  meanX <- mean(X)

  h2n <- compute_h2n(length(X), rho, window_name, window_params, custom_window)
  covVals <- c()

  n <- length(X)
  taperVals_t <- compute_taper(((1:n) - 1/2)/n, rho, window_name, window_params, custom_window)

  for(i in 0:maxLag) {
    taperVals_h <- compute_taper((((1:(n-i)) - 1/2) +i) / n, rho, window_name, window_params, custom_window)

    covVals <- c(covVals, compute_tapered_cov_single(X, meanX, i, h2n, taperVals_t, taperVals_h))
  }

  return(covVals)
}

