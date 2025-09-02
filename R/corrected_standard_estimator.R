#' Kernel Correction of the Standard Estimator.
#'
#' @description
#' This function computes the standard autocovariance estimator and applies kernel correction to it,
#' \deqn{\widehat{C}_{T}^{(a)}(h) = \widehat{C}(h) a_{T}(h), } where \eqn{a_{T}(h) := a(h / N_{T}).}
#' It uses a kernel \eqn{a(\cdot)} which decays or vanishes to zero (depending on the type of kernel) where \eqn{a(0) = 1.}
#' The rate or value at which the kernel vanishes is \eqn{N_{T}}, which is recommended to be of order \eqn{0.1 N}, where \eqn{N} is the length of the observation window, however, one may need to play with this value.
#'
#' @details
#' The aim of this estimator is gradually bring the estimated values to zero through the use of a kernel multiplier. This can be useful when estimating an
#' autocovariance function that is short-range dependent as estimators can have large fluctuations as the lag increases, or to deal with the wave artefacts for large lags, see Bilchouris and Olenko (2025).
#' This estimator can be positive-definite depending on whether the choice of \eqn{\widehat{C}(\cdot)} and \eqn{a} are chosen to be positive-definite or not.
#'
#' @references
#' Yaglom, AM (1987). Correlation Theory of Stationary and Related Random Functions. Volume I: Basic Results. Springer New York. https://doi.org/10.1007/978-1-4612-4628-2
#'
#' Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society (Vol. 54, Issue 1). https://doi.org/10.17713/ajs.v54i1.1975
#'
#' @param X A vector representing observed values of the time series.
#' @param kernel_name The name of the [kernel] function to be used. Possible values are:
#' gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, cauchy.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel] for parameters.
#' In the case of gaussian, wave, rational_quadratic, spherical and circular, \code{N_T} takes the place of \eqn{\theta}.
#' For kernels that require parameters other than \eqn{\theta}, such as the Matern kernel, those parameters are passed.
#' @param N_T The range at which the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a lag small estimation lag.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param custom_kernel If a custom kernel is to be used or not. Defaults to \code{FALSE}. See examples.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the kernel corrected autocovariance estimates.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' corrected_est(X, "gaussian")
#'
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(corrected_est(Y, "bessel_j",
#'      kernel_params=c(0, 1), N_T=0.2*length(Y)))
#'
#' # Custom kernel
#' my_kernel <- function(x, theta, params) {
#'   stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
#'   return(sapply(x, function(t) ifelse(t == 0, 1,
#'          ifelse(t == Inf, 0,
#'          (sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
#' }
#' plot(corrected_est(Y,
#'      my_kernel, kernel_params=c(2, 0.25), custom_kernel = TRUE))
corrected_est <- function(X, kernel_name, kernel_params = c(), N_T = 0.1 * length(X), pd = TRUE, maxLag = length(X) - 1, type = "autocovariance", meanX = mean(X), custom_kernel = FALSE) {
  stopifnot(is.logical(custom_kernel), length(X) > 0, is.vector(X), is.numeric(X), is.numeric(N_T),
            N_T > 0, is.numeric(meanX), is.logical(pd), is.numeric(maxLag), maxLag >= 0,
            maxLag <= (length(X) - 1), maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            type %in% c('autocovariance', 'autocorrelation'))

  retVec <- standard_est(X, pd = pd, maxLag = maxLag, type='autocovariance', meanX = meanX)
  if(!custom_kernel) {
    stopifnot(kernel_name %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))
    retVec <- retVec * sapply(seq(0, maxLag, by=1), function(t) kernel(t, kernel_name, c(N_T, kernel_params)))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    return(retVec)
  }

  if(custom_kernel) {
    stopifnot(exists(quote(kernel_name)))
    retVec <- retVec * sapply(seq(0, maxLag, by=1), function(t) kernel_name(t, N_T, kernel_params))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    return(retVec)
  }

  return("Something went wrong in `corrected_est`.")
}

#' Kernel Correction for an Estimated Autocovariance Function.
#'
#' @description
#' This function applies kernel correction to an estimated autocovariance function,
#' \deqn{\widehat{C}_{T}^{(a)}(h) = \widehat{C}(h) a_{T}(h),} where \eqn{a_{T}(h) := a(h / N_{T}).}
#' It uses a kernel \eqn{a(\cdot)} which decays or vanishes to zero (depending on the type of kernel) where \eqn{a(0) = 1.}
#' The rate or value at which the kernel vanishes is \eqn{N_{T}}, which is recommended to be of order \eqn{0.1 N}, where \eqn{N} is the length of the observation window, however, one may need to play with this value.
#'
#' @param estCov A vector whose values are an estimate autocovariance function.
#' @param kernel_name The name of the [kernel] function to be used. Possible values are:
#' gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, cauchy.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel] for parameters.
#' In the case of gaussian, wave, rational_quadratic, spherical and circular, \code{N_T} takes the place of \eqn{\theta}.
#' For kernels that require parameters other than \eqn{\theta}, such as the Matern kernel, those parameters are passed.
#' @param N_T The range at which the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a lag small estimation lag.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(estCov) - 1}.
#' @param custom_kernel If a custom kernel is to be used or not. Defaults to \code{FALSE}. See the examples of [corrected_est] for usage.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the kernel corrected autocovariance estimates.
#' @export
#'
#' @examples
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' cov_est <- standard_est(Y)
#' plot(cov_est)
#' plot(kernel_est(cov_est,
#'      "bessel_j", kernel_params=c(0, 1), N_T=0.2*length(Y)))
kernel_est <- function(estCov, kernel_name, kernel_params = c(), N_T = 0.1 * length(estCov), maxLag = length(estCov) - 1, type = "autocovariance", custom_kernel = FALSE) {
  stopifnot(is.logical(custom_kernel), length(estCov) > 0, is.vector(estCov), is.numeric(estCov), is.numeric(N_T), N_T > 0,
            is.numeric(maxLag), maxLag >= 0, maxLag <= (length(estCov) - 1), maxLag %% 1 == 0,
            type %in% c('autocovariance', 'autocorrelation'))

  if(!custom_kernel) {
    stopifnot(kernel_name %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))

    estCov <- estCov[1:(maxLag+1)] * sapply(seq(0, maxLag, by=1), function(t) kernel(t, kernel_name, c(N_T, kernel_params)))

    if(type == 'autocorrelation') {
      estCov <- estCov / estCov[1]
    }

    return(estCov)
  }

  if(custom_kernel) {
    stopifnot(exists(quote(kernel_name)))
    estCov <- estCov[1:(maxLag+1)] * sapply(seq(0, maxLag, by=1), function(t) kernel_name(t, N_T, kernel_params))

    if(type == 'autocorrelation') {
      estCov <- estCov / estCov[1]
    }

    return(estCov)
  }

  return("Something went wrong in `kernel_est`.")

}

#' Solve Linear Shrinking
#'
#' This is an objective function used to select \eqn{\lambda \in [0, 1]} in linear shrinking, see [shrinking].
#'
#' @param par The initial parameter  used in the maximisation process.
#' @param corr_mat The autocorrelation matrix of the considered time series.
#' @param target A shrinkage target matrix used in the shrinking process. This should only be used if you wish to use a specific matrix as the target.
#'
#' @references
#' Devlin, S. J., Gnanadesikan R. & Kettenring, J. R. (1975). Robust Estimation and Outlier Detection with Correlation Coefficients. Biometrika, 62(3), 531-545. 10.1093/biomet/62.3.531
#'
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965–984. 10.1080/03610928308831068
#'
#' @return A numeric value that is either equal to \eqn{-}\code{par} or 1.
#' @export
#'
#' @examples
#' estCorr <- c(1, 0.5, 0)
#' corr_mat <- cyclic_matrix(estCorr)
#' solve_shrinking(0.5, corr_mat, diag(length(estCorr)))
solve_shrinking <- function(par, corr_mat, target) {
  stopifnot(is.numeric(par), !is.na(par), length(par) == 1, is.numeric(corr_mat), !any(is.na(corr_mat)), is.numeric(target), !any(is.na(target)),
            all(dim(corr_mat) == dim(target)))
  adj_corr_mat <- (par[1] * corr_mat) + (1 - par[1]) * target
  isPd <- check_pd(adj_corr_mat)
  if(isPd) {
    return(-par[1])
  }
  return(1)
}

#' Linear Shrinking
#'
#' This function corrects an autocovariance/autocorrelation function estimate via linear shrinking of the autocorrelation matrix.
#'
#' @details
#' This function corrects an autocovariance/autocorrelation function estimate via linear shrinking of the autocorrelation matrix.
#' The shrunken autocorrelation matrix is computed as follows
#' \deqn{\widetilde{R} = \lambda R + (1 - \lambda)I_{p},}
#' where \eqn{\widetilde{R}} is the shrunken autocorrelation matrix, \eqn{R} is the original autocorrelation matrix, \eqn{\lambda \in [0, 1],} and \eqn{I_{p}} is the \eqn{p\times p} identity matrix.
#' \eqn{\lambda} is chosen in such a away that largest value which still results in a positive-definite matrix.
#' The shrunken matrix will be positive-definite.
#'
#' @references
#' Devlin, S. J., Gnanadesikan R. & Kettenring, J. R. (1975). Robust Estimation and Outlier Detection with Correlation Coefficients. Biometrika, 62(3), 531-545. 10.1093/biomet/62.3.531
#'
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965–984. 10.1080/03610928308831068
#'
#' @param estCov A vector whose values are an estimate autocovariance/autocorrelation function.
#' @param return_matrix A boolean determining whether the shrunken matrix or the corresponding vector is returned. If \code{FALSE}, it returns a vector whose values are the shrunken autocorrelation function. Defaults to \code{FALSE}.
#' @param target A shrinkage target matrix used in the shrinking process. This should only be used if you wish to use a specific matrix as the target.
#'
#' @return A vector with values of the shrunken autocorrelation function or the corresponding matrix (depending on \code{return_matrix}).
#' @export
#'
#' @examples
#' estCorr <- c(1, 0.8, 0.5, -1.2)
#' shrinking(estCorr)
#' target <- diag(length(estCorr))
#' shrinking(estCorr, TRUE, target)
shrinking <- function(estCov, return_matrix = FALSE, target = NULL) {
  stopifnot(length(estCov) > 0, is.vector(estCov), is.numeric(estCov), !any(is.na(estCov)), is.logical(return_matrix))
  estCorr <- estCov / estCov[1]
  corr_mat <- cyclic_matrix(estCorr)

  if(is.null(target)) {
    target <- diag(nrow(corr_mat))
  }
  else {
    stopifnot(is.matrix(target), is.numeric(target),  !any(is.na(target)))
  }

  optimRes <- optim(0.5, solve_shrinking, corr_mat = corr_mat, target = target, method = 'Brent', lower = 0, upper = 1)
  if(optimRes$convergence != 0) {
    stop(paste0("Failed to converge. Convergence code: ", optimRes$convergence))
  }

  adj_corr_mat <- (optimRes$par * corr_mat) + (1 - optimRes$par) * diag(nrow(corr_mat))

  if(return_matrix) {
    return(list('shrunken'=adj_corr_mat, 'lambda'=optimRes$par))
  }

  return(list('shrunken'=adj_corr_mat[1, ], 'lambda'=optimRes$par))
}
