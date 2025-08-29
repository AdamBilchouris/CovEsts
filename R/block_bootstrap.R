#' Random Block Locations
#'
#' This function performs random sampling to obtian random block locations for block bootstrap.
#'
#' @details
#' This function performs random sampling to obtian random block locations for block bootstrap.
#' If \code{type = 'moving'}, the set \eqn{\{1, \dots, N - \ell  + 1\}} is randomly sampled, with replacement, \eqn{k} times to obtain random block locations for moving block bootstrap.
#' If \code{type = 'circular'}, the set \eqn{\{1, \dots, N\}} is randomly sampled, with replacement, \eqn{k} times to obtain random block locations for moving block bootstrap.
#'
#' @param N The length of the observation window.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @return A vector of length \code{k} whose values are random block locations.
#' @export
#'
#' @examples
#' starting_locs(4, 2, 2)
starting_locs <- function(N, l, k, boot_type = 'moving') {
  stopifnot(is.numeric(N), N >= 1, !is.na(N), is.numeric(l), l >= 1, !is.na(l),
            is.numeric(k), k >= 1, !is.na(k), boot_type %in% c('moving', 'circular'))

  sampleVec <- switch(boot_type,
                      moving = 1:(N - l + 1),
                      circular = 1:N)

  return(sample(sampleVec, k, replace = T))
}

#' Block Bootstrap Samples
#'
#' This function generates block bootstrap samples for either moving block bootstrap or circular bootstrap.
#'
#' @details
#' This function generates a block bootstrap sample for a time series \eqn{X.}
#' For the moving block bootstrap and circular bootstrap procedures see [block_bootstrap] and the included references.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @param X A vector representing observed values of the time series.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#'
#' @return A vector of length \code{length(X)} whose values are a bootstrapped time series.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' bootstrap_samples(X, 2, 3)
bootstrap_samples <- function(X, l, k, boot_type = 'moving') {
  stopifnot(is.numeric(X), length(X) >= 1, !is.na(X), is.numeric(l), l >= 1, !is.na(l),
            is.numeric(k), k >= 1, !is.na(k), boot_type %in% c('moving', 'circular'))
  N <- length(X)
  Y <- c()
  S <- starting_locs(N, l, k, boot_type)

  if(boot_type == 'moving') {
    for(j in 1:k) {
      Y <- c(Y, X[(S[j]):(S[j] + l - 1)])
    }
  }

  else if(boot_type == 'circular') {
    for(j in 1:k) {
      indices <- ( ((S[j]:(S[j] + l - 1)) - 1) %% N) + 1
      Y <- c(Y, X[indices])
    }
  }

  Y <- Y[1:N]
  return(Y)
}

#' Block Bootstrap
#'
#' This function performs block bootstrap (moving or circular) to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#'
#' @details
#' This function performs block bootstrap to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#'
#' Moving block bootstrap can be described as follows.
#' For some times series \eqn{X_{1}, X_{2}, \dots, X_{n},} construct \eqn{k \in Z} overlapping blocks of the form
#' \eqn{B_{i} = (X(i), \dots, X(i + \ell - 1)),} where \eqn{\ell \in [1, n]} is the block length.
#' Randomly sample, with replacement, from the discrete uniform distribution with on \eqn{\{1, \dots, n - \ell + 1\}} to obtain a set of random block locations \eqn{I_{1}, \dots, I_{k}.}
#' Construct a bootstrapped time series \eqn{B_{1}^{\ast}, B_{2}^{\ast}, \dots, B_{k}^{\ast}} where \eqn{B_{i}^{\ast} = B_{I_{i}}.}
#' The bootstrapped time series is truncated to have an observation window of size \eqn{n,} and will be of the form \eqn{X^{\ast}(1), \dots , X^{\ast}(n).}
#' The autocovariance function is then computed for the bootstrapped time series.
#'
#' An alternative to moving block bootstrap is circular block bootstrap.
#' Circular block bootstrap considers the time series to wrap like a circle, that is, the observation at \eqn{n + i \equiv i \; (\bmod \; n).}
#' For example, for the block \eqn{B_{n - \ell + 2}}, we obtain \eqn{(X(n - \ell + 2) , \dots , X(n + 1)) \equiv (X(n - \ell + 2) , \dots , X(1)) \; (\bmod \; n).}
#' When performing random sampling to obtain starting locations, the set \eqn{\{1, \dots, n\}} is now considered.
#' The procedure for constructing the bootstrap time series after constructing the blocks and selecting the starting indices is the same as moving block bootstrap.
#'
#' This process is repeated \code{n_bootstrap} times to obtain many estimates of the autocovariance function using the bootstrapped time series, for which the average autocovariance function
#' and the \eqn{(1 - \alpha)\%} confidence intervals are constructed.
#'
#' The choice of the block length, \eqn{\ell ,} depends on the degree of dependence present in the time series. If the time series has a high degree of dependence, a larger block size should be chosen to ensure the dependency structure is maintained within the block.
#'
#' Any estimator of the autocovariance function can be used in this function, including a custom estimator not in the package, see the examples.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @param X A vector representing observed values of the time series.
#' @param maxLag The maximum lag to compute the estimated autocovariance function at.
#' @param x A vector of indices. Defaults to the sequence \code{1:length(X)}.
#' @param n_bootstrap The number of times to run moving block bootstrap. Defaults to 100.
#' @param l The block length considered for bootstrap. Defaults to \eqn{\lceil N \rceil^{1/3}}, where \eqn{N} is the length of the observation window.
#' @param estimator The function name of the estimator to use. Defaults to \code{standard\_est}.
#' @param alpha The quantile used to compute the \eqn{(1 - \alpha)\%} confidence interval. Defaults to \eqn{0.05.}
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#' @param plot A boolean determining whether a plot should be created. By default, no plot is created.
#' @param boot_mat A boolean determining whether a bootstrap matrix should be returned or not. By default, no matrix is returned.
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
block_bootstrap <- function(X, maxLag, x = 1:length(X), n_bootstrap = 100, l = ceiling(length(X)^(1/3)), estimator = standard_est, alpha = 0.05, boot_type = 'moving', plot = FALSE, boot_mat = FALSE, ...) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(x), length(x) == length(X),
            is.numeric(n_bootstrap), n_bootstrap > 0, is.numeric(l), length(l) == 1, l > 0, l <= length(X),
            l %% 1 == 0, exists(quote(estimator)), is.numeric(alpha), alpha <= 1, alpha >=0,
            boot_type %in% c('moving', 'circular'), is.logical(plot), is.logical(boot_mat))

  N <- length(X)
  k <- ceiling(N / l)
  acf_mat <- matrix(NA, ncol = maxLag + 1, nrow = n_bootstrap)

  for (i in 1:n_bootstrap) {
    Y <- bootstrap_samples(X, l, k, boot_type)
    # If the time series is too small, the same block may be sampled k times.
    while(length(unique(Y)) == 1) {
      Y <- bootstrap_samples(X, l, k, boot_type)
    }
    acf_mat[i, ] <- estimator(Y, maxLag = maxLag, ...)
  }

  lower_vals <- apply(acf_mat, 2, quantile, alpha / 2)
  upper_vals <- apply(acf_mat, 2, quantile, 1 - (alpha / 2))

  original_acf <- estimator(X, ..., maxLag = maxLag)
  average_acf <- colMeans(acf_mat)

  if(plot) {
    plot(x[1:(maxLag + 1)], original_acf, type = "l", ylim = c(-0.6, 1.1), lwd = 2)
    lines(x[1:(maxLag + 1)], average_acf, lty = 2, col = 2, lwd = 2)
    polygon(c(rev(0:maxLag), 0:maxLag), c(rev(upper_vals), lower_vals), col = rgb(0.9, 0.9, 0.9, 0.6))
  }

  if(boot_mat) {
    return(list(average_acf = average_acf, acf_mat = acf_mat, confidence=cbind(lower_vals, upper_vals)))
  }
  return(list(average_acf = average_acf, confidence=cbind(lower_vals, upper_vals)))
}
