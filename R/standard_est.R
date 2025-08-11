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
#' The positive-definite estimator must be used cautiously when estimating over all lags as the sum of all values of the autocorrelation function equals to \eqn{-1/2}.
#' For the nonpositive-definite estimator a similar constant summation property holds.
#'
#' @references Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society 54(1), 112-137. 10.17713/ajs.v54i1.1975
#'
#' @param X A vector representing observed values of the time series.
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
#' standard_est(X, pd = FALSE, maxLag = 2, meanX = mean(X))
standard_est <- function(X, pd = TRUE, maxLag = length(X) - 1, type = "autocovariance", meanX = mean(X)) {
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
#' This function computes an estimated semivariogram using an estimated autocovariance function.
#'
#' @details
#' The semivariogram, \eqn{\gamma(h)} and autocovariance function, \eqn{C(h)}, under the assumption of weak stationarity are related as follows:
#' \deqn{
#' \gamma(h) = C(0) - C(h) .
#' }
#'
#' When an empirical autocovariance function is considered instead, this relation does not necessarily hold, however,
#' it can be used to obtain a function that is close to a semivariogram, see Bilchouris and Olenko (2025).
#'
#' @references Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society 54(1), 112-137. 10.17713/ajs.v54i1.1975
#'
#' @param estCov A vector whose values are an estimate autocovariance function.
#'
#' @return A vector whose values are an estimate of the semivariogram.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' estCov <- standard_est(X, meanX=mean(X), maxLag = 2, pd=FALSE)
#' to_vario(estCov)
to_vario <- function(estCov) {
  stopifnot(length(estCov) > 0, is.vector(estCov), is.numeric(estCov), !any(is.na(estCov)))
  return(estCov[1] - estCov)
}

#' Computes the Standard Estimator of the Autocovariance Function.
#'
#' This function computes the partial autocorrelation function from an estimated autocovariance or autocorrelation function.
#'
#' @details
#' This function is a translation of the 'uni_pacf' function in src/library/stats/src/pacf.c of the R source code which is an implementation of the Durbin–Levinson algorithm.
#'
#' @param estCov A numeric vector representing an estimated autocovariance or autocorrelation function.
#'
#' @return A vector whose values are an estimate partial autocorrelation function.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' to_pacf(standard_est(X, pd = FALSE, maxLag = 2, meanX = mean(X)))
to_pacf <- function(estCov) {
  stopifnot(length(estCov) > 0, is.vector(estCov), is.numeric(estCov))

  # This function is a translation of the 'uni_pacf' function in src/library/stats/src/pacf.c of the R source code.
  n <- length(estCov)
  pacf_vec <- c(estCov[2])
  phi <- c(estCov[2], rep(NA, n - 2))
  temp <- rep(NA, n - 1)

  for(i in 3:(n - 1)) {
    numer <- estCov[i + 1]
    denom <- 1

    for(k in 1:(i - 1)) {
      numer <- numer - (phi[k] * estCov[i - k + 1])
      denom <- denom - (phi[k] * estCov[k + 1])
    }

    pacf_val <- numer / denom
    pacf_vec <- c(pacf_vec, pacf_val)
    phi[i] <- pacf_val

    if (i == (n - 1)) {
      break
    }

    # Create a vector of phi_{n - 1, n - k}
    # This is used in the next calculation.
    for(k in 1:(i - 1)) {
      temp[i - k] <- phi[k]
    }

    # update the phi for next run
    # phi_{n, k} <- phi_{n - 1, k} - (phi_{n, n} * phi_{n - 1, n - k})
    for(k in 1:(i - 1)) {
      phi[k] <- phi[k] - (pacf_val * phi[i - k])
    }
  }

  return(pacf_vec)
}
