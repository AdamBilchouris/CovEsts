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
#' @param window The window function to use, one of
#' "tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch".
#' @param window_params Any paramters of the window function.
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
compute_taper_single <- function(x, rho, window, window_params=c(), custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) == 1, x >= 0, x <= 1, is.numeric(rho), rho > 0)

  if(custom_window) {
    if(x >= 0 && x < ((1/2) * rho)) {
      return(get(window)((2*x) / rho, window_params))
    }
    else if(x >= ((1/2) * rho) && x <= (1/2)) {
      return(1)
    }
    else if(x > (1/2) && x <= 1) {
      x <- 1 - x
      if(x >= 0 && x < ((1/2) * rho)) {
        return(get(window)((2*x) / rho, window_params))
      }
      else if(x >= ((1/2) * rho) && x <= (1/2)) {
        return(1)
      }
    }
  }

  else {
    stopifnot(window %in% c("tukey", "triangular", "power_sine", "blackman_window", "hann_poisson", "welch"))
    if(x >= 0 && x < ((1/2) * rho)) {
      return(get(paste0(window, "_window"))((2*x) / rho, window_params[1]))
    }
    else if(x >= ((1/2) * rho) && x <= (1/2)) {
      return(1)
    }
    else if(x > (1/2) && x <= 1) {
      x <- 1 - x
      if(x >= 0 && x < ((1/2) * rho)) {
        return(get(paste0(window, "_window"))((2*x) / rho, window_params[1]))
      }
      else if(x >= ((1/2) * rho) && x <= (1/2)) {
        return(1)
      }
    }
  }
  return(NA)
}

compute_h2n <- function(n, rho, window, window_params=c(), custom_window = FALSE) {
  stopifnot(is.numeric(n), n >= 1)
  sSeq <- 1:n
  hSeq <- compute_taper(((sSeq - 1/2) / n), rho, window, window_params, custom_window)
  return(sum(hSeq))
}

compute_taper <- function(x, rho, window, window_params, custom_window = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), is.numeric(rho),
            length(rho) == 1, is.logical(custom_window))
  retTaper <- c()
  for(i in 1:length(x)) {
    retTaper <- c(retTaper, compute_taper_single(x[i], rho, window, window_params, custom_window))
  }
  return(retTaper)
}

compute_tapered_cov_single <- function(X, meanX, h, rho, h2n, taperVals_t, taperVals_h, window, window_params = c(), custom_window = FALSE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(meanX), length(meanX) == 1,
            is.numeric(h), length(h) == 1, is.numeric(rho), length(rho) == 1,
            is.numeric(h2n), length(h2n) == 1, is.logical(custom_window))
  covVals <- c()
  n <- length(X)

  for(i in 1:(n - h)) {
    covVals <- c(covVals, (X[i] - meanX) * (X[i + h] - meanX) * taperVals_t[i] * taperVals_h[i])
    # covVals <- c(covVals, (X[i] - meanX) * (X[i + h] - meanX) *
                   # (compute_taper((i - 1/2)/n, rho, window, window_params, custom_window) * compute_taper((i + h - 1/2)/n, rho, window, window_params, custom_window)))
  }
  return(sum(covVals) / h2n)
}

compute_tapered_cov <- function(X, maxLag, rho, window, window_params = c(), custom_window = FALSE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            is.numeric(rho), length(rho) == 1, is.logical(custom_window))

  meanX <- mean(X)

  h2n <- compute_h2n(length(X), rho, window, window_params, custom_window)
  covVals <- c()

  n <- length(X)
  taperVals_t <- compute_taper(((1:n) - 1/2)/n, rho, "tukey", c(), F)

  for(i in 0:maxLag) {
    taperVals_h <- compute_taper((((1:(n-i)) - 1/2) +i) / n, rho, "tukey", c(), F)

    covVals <- c(covVals, compute_tapered_cov_single(X, meanX, i, rho, h2n, taperVals_t, taperVals_h, window, window_params, custom_window))
  }

  return(covVals)
}

