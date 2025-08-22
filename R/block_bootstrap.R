#' Random Block Locations
#'
#' This function randomly samples \eqn{\{1, \dots, N - l  + 1}} \eqn{k} times to obtain block locations for moving block bootstrap.
#'
#' @param N The length of the observation window.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#'
#' @return A vector of length \code{k} whose values are random
#' @export
#'
#' @examples
#' starting_locs(4, 2, 2)
starting_locs <- function(N, l, k) {
  sampleVec <- 1:(N - l + 1)
  return(sample(sampleVec, k, replace = T))
}

#' Title
#'
#' @param X A vector representing observed values of the time series.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#'
#' @return
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' bootstrap(X, 2, 3)
bootstrap <- function(X, l, k) {
  N <- length(X)
  Y <- c()
  S <- starting_locs(N, k)
  for (j in (1:k)) {
    Y <- c(Y, X[(S[j]):(S[j] + l - 1)])
  }
  Y <- Y[1:N]
  return(Y)
}

#' Title
#'
#' @param X A vector representing observed values of the time series.
#' @param maxLag The maximum lag to compute the estimated autocovariance function at.
#' @param n_bootstrap The number of times to run moving block bootstrap. Defaults to 100.
#' @param l The block length considered for bootstrap. Defaults to \eqn{\ceiling{N}^{1/3}}, where \eqn{N} is the length of the observation window.
#' @param estimator The function name of the estimator to use. Defaults to \code{standard\_est}.
#' @param p_quantile The quantile used to compute the \eqn{(1 - \alpha)\%} confidence interval. Defaults to \eqn{0.95.}
#' @param plot A boolean determining whether a plot should be created. By default, no plot is created.
#' @param ... Optional arguments to the chosen estimator. See the examples.
#'
#' @return A list consisting of three items. The f irst is the average estimated autocovariance/autocorrelation function for the bootstrap samples, the second is a matrix of the estimated autocovariance/autocorrelation functions from the bootstrapped samples, and the third is a matrix of confidence intervals for each lag.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' block_bootstrap(X, 4, n_bootstrap = 100, l = 2, type = 'autocorrelation')
#' block_bootstrap(X, 4, n_bootstrap = 100, l = 2, plot =TRUE, type = 'autocovariance')
#' block_bootstrap(X, 4, n_bootstrap = 100, l = 2, estimator=tapered_est, rho = 0.5, window_name = 'blackman', window_params = c(0.16), type='autocorrelation')
block_bootstrap <- function(X, maxLag, n_bootstrap = 100, l = ceiling(length(X)^(1/3)), estimator = standard_est, p_quantile = 0.95, plot = FALSE, ...) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(n_bootstrap), n_bootstrap > 0,
            is.numeric(l), length(l) == 1, l > 0, l <= length(X), l %% 1 == 0, exists(quote(estimator)),
            is.numeric(p_quantile), p_quantile <= 1, p_quantile >=0, is.logical(plot))

  N <- length(X)
  k <- ceiling(N / l)
  acf_mat <- matrix(NA, ncol = maxLag + 1, nrow = n_bootstrap)

  for (i in 1:n_bootstrap) {
    Y <- bootstrap(X, l, k)
    # If the time series is too small, the same block may be sampled k times.
    while(length(unique(Y)) == 1) {
      Y <- bootstrap(X, l, k)
    }
    acf_mat[i, ] <- estimator(Y, maxLag = maxLag, ...)
  }

  lower_quantile <- (1 - p_quantile) / 2
  lower_vals <- apply(acf_mat, 2, quantile, lower_quantile)
  upper_vals <- apply(acf_mat, 2, quantile, p_quantile + lower_quantile)

  original_acf <- estimator(Y, ..., maxLag = maxLag)
  average_acf <- colMeans(acf_mat)

  if(plot) {
    plot(0:maxLag, original_acf, type = "l", ylim = c(-0.6, 1.1), lwd = 2)
    lines(0:maxLag, average_acf, lty = 2, col = 2, lwd = 2)
    polygon(c(rev(0:maxLag), 0:maxLag), c(rev(upper_vals), lower_vals), col = rgb(0.9, 0.9, 0.9, 0.6))
  }

  return(list(average_acf = average_acf, acf_mat = acf_mat, confidence=cbind(lower_vals, upper_vals)))
}
