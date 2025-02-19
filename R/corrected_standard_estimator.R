#' Kernel correction of the standard estimator.
#'
#' This applies kernel correction to the standard estimators.
#' It considers a kernel \eqn{a(\cdot)} which vanishes after a certain distance, \eqn{N_{T}}, which is recommended to be \eqn{0.1 N}, where \eqn{N} is the length of the process.
#' \deqn{\widehat{C}^{(a)}(h) = \widehat{C}(h) a(h / N_{T}).}
#' Also, \eqn{a(0) = 1.}
#'
#' @details
#' The aim of this estimator is gradually bring an estimator to zero through the use of a kernel multiplier. This can be useful when estimating a
#' covariance function that is short-range dependent as estimators can have small fluctuations as the lag increases.
#' This estimator can be positive-definite depending on whether the choice of \eqn{\widehat{C}} and \eqn{a} are chosen to be positive-definite or not.
#'
#' @references
#' Yaglom AM (1987). Correlation Theory of Stationary and Related Random Functions. Volume I: Basic Results. Springer New York. 10.1007/978-1-4612-4628-2.
#'
#' @param X A vector representing the process.
#' @param maxLag The maximum lag to compute the autocovariance function at.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy".
#' @param kernel_params A vector of parameters of the kernel function. See [kernel] for parameters.
#' In the case of "gaussian", "wave", "rational_quadratic", "spherical" and "circular", no parameter is passed as \eqn{\theta} is \code{N_T}.
#' @param N_T The value the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a small range of estimation lags.
#' @param N The length of X.
#' @param meanX The average value of X.
#' @param pd Whether a positive definite estimate should be used.
#' @param type Compute either the 'covariance' or 'correlation'.
#' @param customKernel If a custom kernel is to be used or not. See examples.
#'
#' @return A vector whose values are the estimated autocovariance up to lag maxLag.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_corrected_standard_est(X, 2, "gaussian")
#'
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(compute_corrected_standard_est(Y, length(Y)-1,
#'      "bessel_j", kernel_params=c(0, 1), N_T=0.2*length(Y)))
#'
#' # Custom kernel
#' my_kernel <- function(x, theta, params) {
#'   stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
#'   return(sapply(x, function(t) ifelse(t == 0, 1,
#'          ifelse(t == Inf, 0,
#'          (sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
#' }
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(compute_corrected_standard_est(Y, length(Y)-1,
#'      "my_kernel", kernel_params=c(2, 0.25), customKernel = TRUE))
compute_corrected_standard_est <- function(X, maxLag, kernel_name, kernel_params=c(), N_T=0.1*length(X), N=length(X), meanX=mean(X), pd=TRUE, type='covariance', customKernel = FALSE) {
  stopifnot(is.logical(customKernel), N > 0, length(X) > 0, is.vector(X), is.numeric(X), is.numeric(N_T), N_T > 0, N == length(X), is.numeric(meanX), is.logical(pd),
            is.numeric(maxLag), maxLag >= 0, maxLag <= (N - 1), maxLag %% 1 == 0, type %in% c('covariance', 'correlation'))
  # retVec <- sapply(seq(0, maxLag, by=1), function(tau) standard_est_single(X, tau, N, meanX, pd))
  retVec <- compute_standard_est(X, maxLag, N, meanX, pd, type)
  if(!customKernel) {
    stopifnot(kernel_name %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))

    return(retVec * sapply(seq(0, maxLag, by=1), function(t) get("kernel")(t, kernel_name, c(N_T, kernel_params))))
  }

  if(customKernel) {
    return(retVec * sapply(seq(0, maxLag, by=1), function(t) get(kernel_name)(t, N_T, kernel_params)))
  }

  return("Something went wrong in `compute_corrected_standard_est`.")
}

#' Kernel correction for an estimated autocovariance function.
#'
#' @param cov A vector whose values are an estimate autocovariance function.
#' @param maxLag The maximum lag to compute the autocovariance function at.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy".
#' @param kernel_params A vector of parameters of the kernel function. See [kernel] for parameters.
#' In the case of "gaussian", "wave", "rational_quadratic", "spherical" and "circular", no parameter is passed as \eqn{\theta} is \code{N_T}.
#' @param N_T The value the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a small range of estimation lags.
#' @param customKernel If a custom kernel is to be used or not.  See the examples of [compute_corrected_standard_est] for usage.
#'
#' @return A vector whose values are the estimated autocovariance up to lag maxLag.
#' @export
#'
#' @examples
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' cov_est <- compute_standard_est(Y, length(Y) - 1)
#' plot(compute_kernel_corrected_est(Y, length(Y)-1,
#'      "bessel_j", kernel_params=c(0, 1), N_T=0.2*length(Y)))
compute_kernel_corrected_est <- function(cov, maxLag, kernel_name, kernel_params=c(), N_T=0.1*length(cov), customKernel = FALSE) {
  stopifnot(is.logical(customKernel), length(cov) > 0, is.vector(cov), is.numeric(cov), is.numeric(N_T), N_T > 0,
            is.numeric(maxLag), maxLag >= 0, maxLag <= (length(cov) - 1), maxLag %% 1 == 0)

  if(!customKernel) {
    stopifnot(kernel_name %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))

    return(cov[1:(maxLag+1)] * sapply(seq(0, maxLag, by=1), function(t) get("kernel")(t, kernel_name, c(N_T, kernel_params))))
  }

  if(customKernel) {
    return(cov[1:(maxLag+1)] * sapply(seq(0, maxLag, by=1), function(t) get(kernel_name)(t, N_T, kernel_params)))
  }

  return("Something went wrong in `compute_kernel_corrected_est`.")

}
