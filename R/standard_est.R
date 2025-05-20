#' Computes the autocovariance function using the standard estimator.
#'
#' This function computes the following two estimates of the autocovariance function depending on
#' the parameter \code{pd}.
#'
#' For \code{pd = TRUE}:
#' \deqn{
#' \widehat{C}(\tau) = \frac{1}{N} \sum_{i=1}^{N-\tau} ( X(t_{i}) - \bar{X} ) ( X(t_{i} + \tau) - \bar{X} ) .
#' }
#'
#'For \code{pd = FALSE}:
#'\deqn{
#' \widehat{C}(\tau) = \frac{1}{N - \tau} \sum_{i=1}^{N-\tau} ( X(t_{i}) - \bar{X} ) ( X(t_{i} + \tau) - \bar{X} ) .
#' }
#'
#' This function will generate autocovariance values for lags \eqn{\{0, \dots, \mbox{maxLag}\}.}
#'
#' The positive-definite estimator must be used cautiously when estimating over all lags as the sum of the autocorrelation function converges to \eqn{-1/2}.
#' This does not occur for the nonpositive-definite estimator.
#'
#' @references Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society (Vol. 54, Issue 1). 10.17713/ajs.v54i1.1975
#'
#' @param X A vector representing the observed values of the process.
#' @param maxLag The maximum lag to compute the autocovariance function at.
#' @param N The length of \code{X}. Defaults to \code{length(X)}.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the estimated autocovariance estimates.
#' @export
#'
#' @importFrom stats acf
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_standard_est(X, 2, length(X), mean(X), pd=FALSE)
compute_standard_est <- function(X, maxLag, N=length(X), meanX=mean(X), pd=TRUE, type='autocovariance') {
  stopifnot(N >= 0, length(X) > 0, is.vector(X), is.numeric(X), N == length(X), is.logical(pd), maxLag >= 0, maxLag <= (N - 1),
            maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX), type %in% c('autocovariance', 'autocorrelation'))
  retVec <- as.vector(stats::acf(X, lag.max = maxLag, type = "covariance", plot = FALSE)$acf)
  if(!pd) {
    retVec <- retVec * (length(X) / (length(X) - 0:maxLag))
  }
  if(type == 'autocorrelation') {
    return(retVec / retVec[1])
  }
  return(retVec)
}

#' Autocovariance to Semivariogram
#'
#' This function computes an approximate estimated semivariogram from an an estimated autocovariance function.
#'
#' @details
#' The semivariogram, \eqn{\gamma(h)} and autocovariance function, \eqn{C(h)}, under the assumption of weak stationarity are related as follows:
#' \deqn{
#' \gamma(h) = C(0) - C(h) .
#' }
#'
#' When an empirical autocovariance function is considered instead, this relation does not necessarily hold, however,
#' it can be used to obtain a function that is close to a semivariogram.
#'
#' @param cov A vector whose values are an estimate autocovariance function.
#'
#' @return A vector whose values are an approximate estimate to the semivariogram.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' estCov <- compute_standard_est(X, 2, length(X), mean(X), pd=FALSE)
#' to_vario(estCov)
to_vario <- function(cov) {
  stopifnot(length(cov) > 0, is.vector(cov), is.numeric(cov))
  return(cov[1] - cov)
}
