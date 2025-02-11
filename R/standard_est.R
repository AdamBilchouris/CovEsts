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
#' This function will generate autocovariance values for lags \eqn{\{0, \dots, \mbox{upperTau}\}.}
#'
#' @param X A vector of observations representing a time series.
#' @param upperTau The maximum upper lag to compute the autocovariance function. upperTau \eqn{< N-1.}
#' @param N The length of the vector X.
#' @param meanX The average value of the observations.
#' @param pd Whether a positive definite estimate should be used.
#' @param type Whether the autocovariance or autocorrelation should be computed, options: 'covariance', 'correlation'.
#'
#' @return A vector whose values are the estimated autocovariance up to lag upperTau.
#' @export
#'
#' @importFrom stats acf
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_standard_est(X, 2, length(X), mean(X), pd=FALSE)
compute_standard_est <- function(X, upperTau, N=length(X), meanX=mean(X), pd=TRUE, type='covariance') {
  stopifnot(N >= 0, length(X) > 0, is.vector(X), is.numeric(X), N == length(X), is.logical(pd), upperTau >= 0, upperTau <= (N - 1),
            upperTau %% 1 == 0, type %in% c('covariance', 'correlation'))
  # retVec <- sapply(seq(0, upperTau, by=1), function(tau) standard_est_single(X, tau, N, meanX, pd))
  retVec <- as.vector(stats::acf(X, lag.max = upperTau, type = "covariance", plot = FALSE)$acf)
  if(!pd) {
    retVec <- retVec * (length(X) / (length(X) - 0:upperTau))
  }
  if(type == 'correlation') {
    return(retVec / retVec[1])
  }
  return(retVec)
}
