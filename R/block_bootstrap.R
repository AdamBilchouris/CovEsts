#' Random Block Locations
#'
#' This function randomly samples \eqn{\{1, \dots, N - l  + 1\}} \eqn{k} times to obtain random block locations for moving block bootstrap.
#'
#' @param N The length of the observation window.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#'
#' @references
#' Chapter 2.5 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#'
#' @return A vector of length \code{k} whose values are random block locations.
#' @export
#'
#' @examples
#' starting_locs(4, 2, 2)
starting_locs <- function(N, l, k) {
  sampleVec <- 1:(N - l + 1)
  return(sample(sampleVec, k, replace = T))
}

#' Moving Block Bootstrap Samples
#'
#' This function generates moving block bootstrap samples.
#'
#' @details
#' This function generates a moving block bootstrap sample for a time series \eqn{X.}
#' For some times series \eqn{X_{1}, X_{2}, \dots, X_{n},} construct \eqn{k \in Z} overlapping blocks of the form
#' \eqn{B_{i} = (X_{i}, \dots, X_{i + \ell - 1}),} where \eqn{\ell \in [1, n]} is the block length.
#' Randomly sample, with replacement, from the discrete uniform distribution with on \eqn{\{1, \dots, n - \ell + 1\}} to obtain a set of random block locations \eqn{I_{1}, \dots, I_{k}.}
#' Construct a bootstrapped time series \eqn{X_{1}^{\ast}, X_{2}^{\ast}, \dots, X_{n}^{\ast}} where \eqn{X_{i}^{\ast} = B_{I_{i}}.}
#' The bootstrapped time series is truncated to have an observation window of size \eqn{n.}
#'
#' @references
#' Chapter 2.5 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' @param X A vector representing observed values of the time series.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#'
#' @return A vector of length \code{length(X)} whose values are a bootstrapped time series.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' bootstrap_samples(X, 2, 3)
bootstrap_samples <- function(X, l, k) {
  N <- length(X)
  Y <- c()
  S <- starting_locs(N, k)
  for(j in 1:k) {
    Y <- c(Y, X[(S[j]):(S[j] + l - 1)])
  }
  Y <- Y[1:N]
  return(Y)
}

#' Moving Block Bootstrap
#'
#' This function performs moving block bootstrap to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#'
#' @details
#' This function performs moving block bootstrap to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#' Moving block bootstrap can be described as follows.
#' For some times series \eqn{X_{1}, X_{2}, \dots, X_{n},} construct \eqn{k \in Z} overlapping blocks of the form
#' \eqn{B_{i} = (X_{i}, \dots, X_{i + \ell - 1}),} where \eqn{\ell \in [1, n]} is the block length.
#' Randomly sample, with replacement, from the discrete uniform distribution with on \eqn{\{1, \dots, n - \ell + 1\}} to obtain a set of random block locations \eqn{I_{1}, \dots, I_{k}.}
#' Construct a bootstrapped time series \eqn{X_{1}^{\ast}, X_{2}^{\ast}, \dots, X_{n}^{\ast}} where \eqn{X_{i}^{\ast} = B_{I_{i}}.}
#' The bootstrapped time series is truncated to have an observation window of size \eqn{n.}
#' The autocovariance function is then computed for the bootstrapped time series.
#'
#' This process is repeated \code{n\_bootstrap} times to obtain many estimates of the autocovariance function using the bootstrapped time series, for which the average autocovariance function
#' and the \eqn{(1 - \alpha)\%} confidence intervals are constructed.
#'
#' Any estimator of the autocovariance function can be used in this function, including a custom estimator not in the package, see the examples.
#'
#' @references
#' Chapter 2.5 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' @param X A vector representing observed values of the time series.
#' @param maxLag The maximum lag to compute the estimated autocovariance function at.
#' @param n_bootstrap The number of times to run moving block bootstrap. Defaults to 100.
#' @param l The block length considered for bootstrap. Defaults to \eqn{\lceil N \rceil^{1/3}}, where \eqn{N} is the length of the observation window.
#' @param estimator The function name of the estimator to use. Defaults to \code{standard\_est}.
#' @param alpha The quantile used to compute the \eqn{(1 - \alpha)\%} confidence interval. Defaults to \eqn{0.05.}
#' @param plot A boolean determining whether a plot should be created. By default, no plot is created.
#' @param ... Optional named arguments to the chosen estimator. See the examples.
#'
#' @return A list consisting of three items. The first is the average estimated autocovariance/autocorrelation function for the bootstrap samples, the second is a matrix of the estimated autocovariance/autocorrelation functions from the bootstrapped samples, and the third is a matrix of confidence intervals for each lag.
#' @export
#'
#' @importFrom stats quantile
#' @importFrom grDevices rgb
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, type = 'autocorrelation')
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, plot =TRUE, type = 'autocovariance')
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, estimator=tapered_est,
#'     rho = 0.5, window_name = 'blackman', window_params = c(0.16),
#'     type='autocorrelation')
#
#' my_cov_est <- function(X, maxLag) {
#'   n <- length(X)
#'   covVals <- c()
#'   for(h in 0:maxLag) {
#'     covVals_t <- (X[1:(n-h)] - mean(X)) * (X[(1+h):n] - mean(X))
#'     covVals <- c(covVals, sum(covVals_t) / (n - h))
#'   }
#'   return(covVals)
#' }
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, estimator=my_cov_est)
block_bootstrap <- function(X, maxLag, n_bootstrap = 100, l = ceiling(length(X)^(1/3)), estimator = standard_est, alpha = 0.05, plot = FALSE, ...) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(n_bootstrap), n_bootstrap > 0,
            is.numeric(l), length(l) == 1, l > 0, l <= length(X), l %% 1 == 0, exists(quote(estimator)),
            is.numeric(alpha), alpha <= 1, alpha >=0, is.logical(plot))

  N <- length(X)
  k <- ceiling(N / l)
  acf_mat <- matrix(NA, ncol = maxLag + 1, nrow = n_bootstrap)

  for (i in 1:n_bootstrap) {
    Y <- bootstrap_samples(X, l, k)
    # If the time series is too small, the same block may be sampled k times.
    while(length(unique(Y)) == 1) {
      Y <- bootstrap_samples(X, l, k)
    }
    acf_mat[i, ] <- estimator(Y, maxLag = maxLag, ...)
  }

  lower_vals <- apply(acf_mat, 2, quantile, alpha / 2)
  upper_vals <- apply(acf_mat, 2, quantile, 1 - (alpha / 2))

  original_acf <- estimator(X, ..., maxLag = maxLag)
  average_acf <- colMeans(acf_mat)

  if(plot) {
    plot(0:maxLag, original_acf, type = "l", ylim = c(-0.6, 1.1), lwd = 2)
    lines(0:maxLag, average_acf, lty = 2, col = 2, lwd = 2)
    polygon(c(rev(0:maxLag), 0:maxLag), c(rev(upper_vals), lower_vals), col = rgb(0.9, 0.9, 0.9, 0.6))
  }

  return(list(average_acf = average_acf, acf_mat = acf_mat, confidence=cbind(lower_vals, upper_vals)))
}
