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
#' @param kernel_name The name of the [kernel_ec] function to be used. Possible values are:
#' gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, cauchy.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel_ec] for parameters.
#' In the case of gaussian, wave, rational_quadratic, spherical and circular, \code{N_T} takes the place of \eqn{\theta}.
#' For kernels that require parameters other than \eqn{\theta}, such as the Matern kernel, those parameters are passed.
#' @param N_T The range at which the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a lag small estimation lag.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param x A vector of lag indices. Defaults to the sequence \code{0:length(X)}. Must be at least as large as \code{maxLag + 1}.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param custom_kernel If a custom kernel is to be used or not. Defaults to \code{FALSE}. See examples.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the autocovariance/autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the \code{type} parameter.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'corrected_est'.}
#' }
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
corrected_est <- function(X, kernel_name = c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"), kernel_params = c(), N_T = 0.1 * length(X), pd = TRUE, maxLag = length(X) - 1, x = 0:length(X), type = c("autocovariance", "autocorrelation"), meanX = mean(X), custom_kernel = FALSE) {
  stopifnot(is.logical(custom_kernel), length(X) > 0, all(!is.na(X)), is.numeric(X), is.numeric(N_T), !is.na(N_T),
            N_T > 0, is.numeric(meanX), !is.na(meanX), is.logical(pd), is.numeric(maxLag), maxLag >= 0, !is.na(maxLag),
            maxLag <= (length(X) - 1), maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            is.numeric(x), all(!is.na(x)), length(x) > 0)

  type <- match.arg(type)

  retVec <- standard_est(X, pd = pd, maxLag = maxLag, type='autocovariance', meanX = meanX)$acf
  if(!custom_kernel) {
    kernel_name <- match.arg(kernel_name)
    retVec <- retVec * kernel_ec(0:maxLag, kernel_name, c(N_T, kernel_params))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    res <- list(acf = retVec, lags = x[1:(maxLag + 1)], est_type = type, est_used = 'corrected_est')
    return(structure(res, class = "CovEsts"))
  }

  if(custom_kernel) {
    stopifnot(is.function(kernel_name))
    # Keep this as an apply statement as the user may not implement a kernel that allows for a vector paramter.
    retVec <- retVec * sapply(seq(0, maxLag, by=1), function(t) kernel_name(t, N_T, kernel_params))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    res <- list(acf = retVec, lags = x[1:(maxLag + 1)], est_type = type, est_used = 'corrected_est')
    return(structure(res, class = "CovEsts"))
  }
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
#' @param kernel_name The name of the [kernel_ec] function to be used. Possible values are:
#' gaussian, exponential, wave, rational_quadratic, spherical, circular, bessel_j, matern, cauchy.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel_ec] for parameters.
#' In the case of gaussian, wave, rational_quadratic, spherical and circular, \code{N_T} takes the place of \eqn{\theta}.
#' For kernels that require parameters other than \eqn{\theta}, such as the Matern kernel, those parameters are passed.
#' @param N_T The range at which the kernel function vanishes at. Recommended to be \eqn{0.1 N} when considering all lags. This parameter may be large for a lag small estimation lag.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(estCov) - 1}.
#' @param x A vector of lag indices. Defaults to the sequence \code{0:length(X)}. Must be at least as large as \code{maxLag + 1}.
#' @param custom_kernel If a custom kernel is to be used or not. Defaults to \code{FALSE}. See the examples of [corrected_est] for usage.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the kernel corrected autocovariance estimates or \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the autocovariance/autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on, inherited from the argument \code{estCov}.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the argument \code{type}.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'kernel_est'.}
#' }
#' If a numeric vector is given for the argument \code{estCov}, then a numeric vector output is given, and if a \code{CovEsts} S3 object is given, a \code{CovEsts} object is given as output.
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
kernel_est <- function(estCov, kernel_name = c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"), kernel_params = c(), N_T = 0.1 * length(estCov), maxLag = length(estCov) - 1, x = 0:length(estCov), type = c("autocovariance", "autocorrelation"), custom_kernel = FALSE) {
  UseMethod("kernel_est")
}

#' @describeIn kernel_est Method for `CovEsts` objects.
#' @export
kernel_est.CovEsts <- function(estCov, kernel_name = c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"), kernel_params = c(), N_T = 0.1 * length(estCov$acf), maxLag = length(estCov$acf) - 1, x = estCov$lags, type = c("autocovariance", "autocorrelation"), custom_kernel = FALSE) {
  est <- kernel_est.default(estCov$acf, kernel_name, kernel_params = kernel_params, N_T = N_T, maxLag = length(estCov$lags) - 1, type = type, custom_kernel = custom_kernel)

  res <- list(acf = est, lags = estCov$lags, est_type = type, est_used = 'kernel_est')
  return(structure(res, class = "CovEsts"))
}

#' @describeIn kernel_est Method for numeric vectors.
#' @export
kernel_est.default <- function(estCov, kernel_name = c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"), kernel_params = c(), N_T = 0.1 * length(estCov), maxLag = length(estCov) - 1, x = 0:length(estCov), type = c("autocovariance", "autocorrelation"), custom_kernel = FALSE) {
  stopifnot(is.logical(custom_kernel), length(estCov) > 0, is.numeric(estCov), !any(is.na(estCov)), is.numeric(N_T), N_T > 0,
            !is.na(N_T), is.numeric(maxLag), maxLag >= 0, maxLag <= (length(estCov) - 1), maxLag %% 1 == 0, !is.na(maxLag))

  type <- match.arg(type)

  if(!custom_kernel) {
    kernel_name <- match.arg(kernel_name)

    estCov <- estCov[1:(maxLag+1)] * kernel_ec(0:maxLag, kernel_name, c(N_T, kernel_params))

    if(type == 'autocorrelation') {
      estCov <- estCov / estCov[1]
    }

    return(estCov)
  }

  if(custom_kernel) {
    stopifnot(is.function(kernel_name))
    # Keep this as an apply statement as the user may not implement a kernel that allows for a vector paramter.
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
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965-984. 10.1080/03610928308831068
#'
#' @return A numeric value that is either equal to \eqn{-}\code{par} or 1.
#'
#' @examples
#' \dontrun{
#' estCorr <- c(1, 0.5, 0)
#' corr_mat <- cyclic_matrix(estCorr)
#' solve_shrinking(0.5, corr_mat, diag(length(estCorr)))
#' }
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
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965-984. 10.1080/03610928308831068
#'
#' @param estCov A vector whose values are an estimate autocovariance/autocorrelation function.
#' @param return_matrix A boolean determining whether the shrunken matrix or the corresponding vector is returned. If \code{FALSE}, it returns a vector whose values are the shrunken autocorrelation function. Defaults to \code{FALSE}.
#' @param target A shrinkage target matrix used in the shrinking process. This should only be used if you wish to use a specific matrix as the target.
#'
#' @return A vector with values of the shrunken autocorrelation function, the corresponding matrix (depending on \code{return_matrix}), or \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the shrunken autocovariance/autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on, inherited from the argument \code{estCov}.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the argument \code{type}.}
#'  \item{\code{est_used}}{The estimator function used, in this case, inherited from the argument \code{estCov}.}
#'  \item{\code{correction_method}}{This value is 'shrinking'.}
#'  \item{\code{lambda}}{The \eqn{\lambda} value obtained during the shrinking process.}
#' }
#' If a numeric vector is given for the argument \code{estCov}, then a numeric vector output is given, and if a \code{CovEsts} S3 object is given, a \code{CovEsts} object is given as output.
#'
#' @export
#'
#' @examples
#' estCorr <- c(1, 0.8, 0.5, -1.2)
#' shrinking(estCorr)
#' target <- diag(length(estCorr))
#' shrinking(estCorr, TRUE, target)
shrinking <- function(estCov, return_matrix = FALSE, target = NULL) {
  UseMethod("shrinking")
}

#' @describeIn shrinking Method for `CovEsts` objects.
#' @export
shrinking.CovEsts <- function(estCov, return_matrix = FALSE, target = NULL) {
  shrunk <- shrinking.default(as.numeric(estCov), return_matrix = return_matrix, target = target)

  res <- list(acf=shrunk[['shrunken']], lags = estCov$lags, est_type = estCov$est_type,
              est_used = estCov$est_used, correction_method = 'shrinking',
              lambda = shrunk[['lambda']])
  return(structure(res, class = "CovEsts"))
}

#' @describeIn shrinking Method for numeric vectors
#' @export
shrinking.default <- function(estCov, return_matrix = FALSE, target = NULL) {
  stopifnot(length(estCov) > 0, is.numeric(estCov), !any(is.na(estCov)), is.logical(return_matrix))
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
