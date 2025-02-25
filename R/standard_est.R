#' Computes the autocovariance function using the standard estimator.
#'
#' This computes the following two variants of the autocovariance function depending on
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
#' @param X A vector representing the process.
#' @param maxLag The maximum lag to compute the autocovariance function at.
#' @param N The length of X.
#' @param meanX The average value of X.
#' @param pd Whether a positive-definite estimate should be used.
#' @param type Compute either the 'covariance' or 'correlation'.
#'
#' @return A vector whose values are the estimated autocovariance.
#' @export
#'
#' @importFrom stats acf
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_standard_est(X, 2, length(X), mean(X), pd=FALSE)
compute_standard_est <- function(X, maxLag, N=length(X), meanX=mean(X), pd=TRUE, type='covariance') {
  stopifnot(N >= 0, length(X) > 0, is.vector(X), is.numeric(X), N == length(X), is.logical(pd), maxLag >= 0, maxLag <= (N - 1),
            maxLag %% 1 == 0, type %in% c('covariance', 'correlation'))
  # retVec <- sapply(seq(0, maxLag, by=1), function(tau) standard_est_single(X, tau, N, meanX, pd))
  retVec <- as.vector(stats::acf(X, lag.max = maxLag, type = "covariance", plot = FALSE)$acf)
  if(!pd) {
    retVec <- retVec * (length(X) / (length(X) - 0:maxLag))
  }
  if(type == 'correlation') {
    return(retVec / retVec[1])
  }
  return(retVec)
}
