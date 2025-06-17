#' Computes the Standard Estimator of the Autocovariance Function.
#'
#' This function computes the following two estimates of the autocovariance function depending on
#' the parameter \code{pd}.
#'
#' For \code{pd = TRUE}:
#' \deqn{
#' \widehat{C}(h) = \frac{1}{N} \sum_{j=1}^{N-h} ( X(j) - \bar{X} ) ( X(j + h) - \bar{X} ) .
#' }
#'
#'For \code{pd = FALSE}:
#'\deqn{
#' \widehat{C}(h) = \frac{1}{N - h} \sum_{j=1}^{N-h} ( X(j) - \bar{X} ) ( X(j + h) - \bar{X} ) .
#' }
#'
#' This function will generate autocovariance values for lags \eqn{h} from the set \eqn{\{0, \dots, \mbox{maxLag}\}.}
#'
#' The positive-definite estimator must be used cautiously when estimating over all lags as the sum of the autocorrelation function converges to \eqn{-1/2}.
#' This does not occur for the nonpositive-definite estimator.
#'
#' @references Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society (Vol. 54, Issue 1). 10.17713/ajs.v54i1.1975
#'
#' @param X A vector with observed values.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the autocovariance estimates.
#' @export
#'
#' @importFrom stats acf
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_standard_est(X, meanX=mean(X), maxLag=2, pd=FALSE)
compute_standard_est <- function(X, meanX=mean(X), maxLag=length(X) - 1, pd=TRUE, type='autocovariance') {
  stopifnot(length(X) > 0, is.vector(X), is.numeric(X), is.logical(pd), maxLag >= 0, maxLag <= (length(X) - 1),
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
#' estCov <- compute_standard_est(X, meanX=mean(X), maxLag = 2, pd=FALSE)
#' to_vario(estCov)
to_vario <- function(cov) {
  stopifnot(length(cov) > 0, is.vector(cov), is.numeric(cov), !any(is.na(cov)))
  return(cov[1] - cov)
}
