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
#' Yaglom, AM (1987). Correlation Theory of Stationary and Related Random Functions. Volume I: Basic Results. Springer New York. 10.1007/978-1-4612-4628-2.
#'
#' Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society (Vol. 54, Issue 1). 10.17713/ajs.v54i1.1975
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
#' compute_corrected_standard_est(X, "gaussian")
#'
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' plot(Y)
#' plot(compute_corrected_standard_est(Y, "bessel_j",
#'      kernel_params=c(0, 1), N_T=0.2*length(Y)))
#'
#' # Custom kernel
#' my_kernel <- function(x, theta, params) {
#'   stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
#'   return(sapply(x, function(t) ifelse(t == 0, 1,
#'          ifelse(t == Inf, 0,
#'          (sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
#' }
#' plot(compute_corrected_standard_est(Y,
#'      my_kernel, kernel_params=c(2, 0.25), custom_kernel = TRUE))
compute_corrected_standard_est <- function(X, kernel_name, kernel_params = c(), N_T = 0.1 * length(X), pd = TRUE, maxLag = length(X) - 1, type = "autocovariance", meanX = mean(X), custom_kernel = FALSE) {
  stopifnot(is.logical(custom_kernel), length(X) > 0, is.vector(X), is.numeric(X), is.numeric(N_T),
            N_T > 0, is.numeric(meanX), is.logical(pd), is.numeric(maxLag), maxLag >= 0,
            maxLag <= (length(X) - 1), maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            type %in% c('autocovariance', 'autocorrelation'))

  retVec <- compute_standard_est(X, pd = pd, maxLag = maxLag, type='autocovariance', meanX = meanX)
  if(!custom_kernel) {
    stopifnot(kernel_name %in% c("gaussian", "exponential", "wave", "rational_quadratic", "spherical", "circular", "bessel_j", "matern", "cauchy"))
    retVec <- retVec * sapply(seq(0, maxLag, by=1), function(t) kernel(t, kernel_name, c(N_T, kernel_params)))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    return(retVec)
  }

  if(custom_kernel) {
    retVec <- retVec * sapply(seq(0, maxLag, by=1), function(t) kernel_name(t, N_T, kernel_params))

    if(type == 'autocorrelation') {
      retVec <- retVec / retVec[1]
    }

    return(retVec)
  }

  return("Something went wrong in `compute_corrected_standard_est`.")
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
#' @param custom_kernel If a custom kernel is to be used or not. Defaults to \code{FALSE}. See the examples of [compute_corrected_standard_est] for usage.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the kernel corrected autocovariance estimates.
#' @export
#'
#' @examples
#' X <- rnorm(1000)
#' Y <- c(X[1], X[2])
#' for(i in 3:length(X)) { Y[i] <- X[i] - 0.3*X[i - 1] - 0.6*X[i - 2] }
#' cov_est <- compute_standard_est(Y)
#' plot(cov_est)
#' plot(compute_kernel_corrected_est(cov_est,
#'      "bessel_j", kernel_params=c(0, 1), N_T=0.2*length(Y)))
compute_kernel_corrected_est <- function(estCov, kernel_name, kernel_params = c(), N_T = 0.1 * length(estCov), maxLag = length(estCov) - 1, type = "autocovariance", custom_kernel = FALSE) {
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
    estCov <- estCov[1:(maxLag+1)] * sapply(seq(0, maxLag, by=1), function(t) kernel_name(t, N_T, kernel_params))

    if(type == 'autocorrelation') {
      estCov <- estCov / estCov[1]
    }

    return(estCov)
  }

  return("Something went wrong in `compute_kernel_corrected_est`.")

}

#' Solve Linear Shrinking
#'
#' This is an objective function used to select \eqn{\lambda} in linear shrinking, see [compute_linear_shrinking].
#'
#' @param par The initial parameter  used in the maximisation process.
#' @param corr_mat The autocorrelation matrix of the considered time series.
#'
#' @references
#' Devlin, S. J., Gnanadesikan R. & Kettenring, J. R. (1975). Robust Estimation ofand Outlier Detection with Correlation Coefficients. Biometrika, 62(3), 531-545. 10.1093/biomet/62.3.531
#'
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965–984. 10.1080/03610928308831068
#'
#' @return A numeric value that is either equal to \eqn{-}\code{par} or 1.
#' @export
#'
#' @examples
#' estCorr <- c(1, 0.5, 0)
#' corr_mat <- create_cyclic_matrix(estCorr)
#' solve_linear_shrinking(0.5, corr_mat)
solve_linear_shrinking <- function(par, corr_mat) {
  adj_corr_mat <- (par[1] * corr_mat) + (1 - par[1]) * diag(nrow(corr_mat))
  isPd <- check_pd(adj_corr_mat[1, ])
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
#' This function corrects an autocovariance/autocorrelation function estimate via linear shrinking of the autocorrelation matrix (see a b).
#' The shrunken autocorrelation matrix is computed as follows
#' \deqn{\widetilde{R} \coloneqq \lambda R + (1 - \lambda)I_{p},}
#' where \eqn{\widetilde{R}} is the shrunken autocorrelation matrix, \eqn{R} is the original autocorrelation matrix, \eqn{\lambda \in [0, 1],} and \eqn{I_{p}} is the \eqn{p\times p} identity matrix.
#' \eqn{\lambda} is chosen in such a away that largest value which still results in a positive-definite matrix.
#' The shrunken matrix will be positive-definite.
#'
#' @references
#' Devlin, S. J., Gnanadesikan R. & Kettenring, J. R. (1975). Robust Estimation ofand Outlier Detection with Correlation Coefficients. Biometrika, 62(3), 531-545. 10.1093/biomet/62.3.531
#'
#' Rousseeuw, P. J. & Molenberghs, G. (1993). Transformation of Non Positive Semidefinite Correlation Matrices. Communications in Statistics - Theory and Methods, 22(4), 965–984. 10.1080/03610928308831068
#'
#' @param estCov A vector whose values are an estimate autocovariance/autocorrelation function.
#' @param return_matrix A boolean determining whether or not the shrunken matrix is given or not. If \code{FALSE}, it returns the shrunken autocorrelation function. Defaults to \code{FALSE}.
#'
#' @return The shrunken autocorrelation function or matrix (depending on \code{return_matrix}).
#' @export
#'
#' @examples
#' estCorr <- c(1, 0.8, 0.5, -1.2)
#' compute_linear_shrinking(estCorr)
#' compute_linear_shrinking(estCorr, TRUE)
compute_linear_shrinking <- function(estCov, return_matrix = FALSE) {
  estCorr <- estCov / estCov[1]
  corr_mat <- create_cyclic_matrix(estCorr)
  optimRes <- optim(0.5, solve_linear_shrinking, corr_mat = corr_mat, method = 'Brent', lower = 0, upper = 1)
  if(optimRes$convergence != 0) {
    stop(paste0("Failed to converge. Convergence code: ", optimRes$convergence))
  }

  adj_corr_mat <- (optimRes$par * corr_mat) + (1 - optimRes$par) * diag(nrow(corr_mat))

  if(return_matrix) {
    return(list('shrunken'=adj_corr_mat, 'lambda'=optimRes$par))
  }

  return(list('shrunken'=adj_corr_mat[1, ], 'lambda'=optimRes$par))
}

#' Compute the Long-Run Covariance Sandwich Estimator
#'
#' @description
#' This function computes the long-run covariance sandwich estimator of a time series, proposed in Shang (2024).
#'
#' @details
#' This function computes the long-run covariance sandwich estimator of a time series.
#' The long-run covariance is a measure of total dependence in a time series as opposed to at a single lag, like the autocovariance.
#'
#' Shang (2024) provides the following estimator (adapted for the nonfunctional time series case)
#' \deqn{
#' \widehat{C}_{b} \coloneqq \sum_{h=-(N-1)}^{N-1} w(h / b) \widehat{C}(h) ,
#' }
#' where \eqn{b > 0} is the bandwith parameter, \eqn{w(\cdot)} is a symmetric window function such that \eqn{w(0) = 1, w(\cdot) \leq 1} and after some value \eqn{\iota} (we selected 1) \eqn{w(\cdot) = 0}.
#' In the paper, \eqn{\widehat{C}(\cdot)} is the standard positive-definite autocovariance estimator (see [compute_standard_est]).
#'
#'
#' @references Shang, H. L. (2024). Bootstrapping Long-Run Covariance of Stationary Functional Time Series. Forecasting 6(1), 138-151. 10.3390/forecast6010008
#'
#' @param X A vector representing observed values of the time series.
#' @param window_name The name of the [window] function to be used. Possible values are:
#' tukey, triangular, power_sine, blackman_window, hann_poisson, welch. Alternatively, a custom window function can be provided, see the example
#' @param window_params A vector of parameters of the window function.
#' @param b Bandwidth parameter, greater than 0. Defaults to \eqn{\lfloor} length(X) \eqn{\rfloor^{1/5}.}
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated long-run covariance. Defaults to \code{length(X) - 1}.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param custom_window If a custom window is to be used or not. Defaults to \code{FALSE}.
#'
#' @return The long-run covariance of the time series.
#' @export
#'
#' @examples
#' X <- 0:4
#' compute_sandwich_est(X)
#' my_window <- function(x, ...) {
#'   retVal <- 1 - x
#'   retVal[abs(x) > 1] <- 0
#'   return(retVal)
#' }
#' compute_sandwich_est(X, window_name=my_window, custom_window = TRUE)
compute_sandwich_est <- function(X, window_name = 'tukey', window_params = c(1), b = floor(length(X)^(1/5)), pd = TRUE, maxLag = length(X) - 1, meanX = mean(X), custom_window = FALSE) {
  stopifnot(length(X) > 0, is.vector(X), is.numeric(X), is.logical(pd), maxLag >= 0, maxLag <= (length(X) - 1),
            maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            length(b) == 1, is.numeric(b), !is.na(b), is.logical(custom_window))
  base_est <- compute_standard_est(X, pd = pd, maxLag = maxLag, type = 'autocovariance', meanX = meanX)
  base_est <- c(rev(base_est[-1]), base_est)
  window_vals <- (0:maxLag) / b
  window_vals <- c(rev(window_vals[-1]), window_vals)

  if(custom_window) {
    return(sum(base_est * window_name(window_vals, window_params)))
  }
  else {
    stopifnot(window_name %in% c("tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"))
    return(sum(base_est * window_symm(window_vals, window_name, window_params)))
  }

}
