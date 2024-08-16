#' Kernel correction of the standard estimator.
#'
#' This applies kernel correction to the standard estimators.
#' It considers a kernel \eqn{a(\cdot)} which vanishes after a certain distancem, \eqn{N_{T}}, which is recommended to be \eqn{0.1 N}, where \eqn{N} is the length of the process.
#' \deqn{\widehat{C}^{(a)}(h) = \widehat{C}(h) a(h / N_{T}).}
#' Also, \eqn{a(0) = 1.}
#'
#' @references
#' Yaglom AM (1987). Correlation Theory of Stationary and Related Random Functions. Volume I: Basic Results. Springer New York. 10.1007/978-1-4612-4628-2.
#'
#' @param X A vector of observations representing a time series.
#' @param upperTau The maximum upper lag to compute the autocovariance function. upperTau \eqn{< N-1.}
#' @param kernel The choice of kernel. Possible values are:
#' "gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy".
#' @param kernel_params If the kernel has any parameters, pass them as a vector.
#' @param N_T The value the kernel function vanishes at. Recommended to be \eqn{0.1 N.}
#' @param N The length of the vector X.
#' @param meanX The average value of the observations.
#' @param pd Whether a positive definite estimate should be used.
#' @param type Whether the autocovariance or autocorrelation should be computed, options: 'covariance', 'correlation'.
#' @param customKernel Whether or not you want to use a custom kernel function. This will allow you to define your own kernel function. See the examples for usage.
#'
#' @return The autocovariance up to lag upperTau for the time series X.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' corrected_standard_estimator(X, 2, "gaussian")
#'
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(corrected_standard_estimator(Y, length(Y)-1,
#'      "bessel_j", kernel_params=c(0, 1), N_T=0.2*length(Y)))
#'
#' # Custom kernel
#' my_kernel <- function(x, theta, params) {
#'   stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
#'   return(sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0, (sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
#' }
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(corrected_standard_estimator(Y, length(Y)-1, "my_kernel", kernel_params=c(2, 0.25), customKernel = TRUE))
corrected_standard_estimator <- function(X, upperTau, kernel, kernel_params=c(), N_T=0.1*length(X), N=length(X), meanX=mean(X), pd=TRUE, type='covariance', customKernel = FALSE) {
  stopifnot(is.logical(customKernel), N >= 0, length(X) > 0, is.vector(X), N == length(X), is.logical(pd), upperTau >= 0, upperTau <= (N - 1),
            type %in% c('covariance', 'correlation'))
  retVec <- sapply(seq(0, upperTau, by=1), function(tau) standard_est_single(X, tau, N, meanX, pd))
  if(!customKernel) {
    stopifnot(kernel %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))

    # Check length of kernel parameters
    # Length 0, e.g. Gaussian, exponential.
    if(length(kernel_params) == 0) {
      return(retVec * sapply(seq(0, upperTau, by=1), function(t) get(paste0("kernel_", kernel))(t, N_T)))
    }

    # Length 1, e.g. Matern
    else if(length(kernel_params) == 1) {
      return(retVec * sapply(seq(0, upperTau, by=1), function(t) get(paste0("kernel_", kernel))(t, N_T, kernel_params[1])))
    }

    # Length 2, e.g. Bessel J, Cauchy
    else {
      return(retVec * sapply(seq(0, upperTau, by=1), function(t) get(paste0("kernel_", kernel))(t, N_T, kernel_params[1], kernel_params[2])))
    }
  }

  if(customKernel) {
    return(retVec * sapply(seq(0, upperTau, by=1), function(t) get(kernel)(t, N_T, kernel_params)))
  }

  return("Something went wrong in `corrected_standard_estimator`.")
}
