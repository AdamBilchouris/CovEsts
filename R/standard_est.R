#' Computes the Standard Estimator of the Autocovariance Function.
#'
#' This function computes the following two estimates of the autocovariance function depending on
#' the parameter \code{pd}, see the Details section.
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
#' @references
#' Bilchouris, A. & Olenko, A (2025). On Nonparametric Estimation of Covariogram. Austrian Statistical Society 54(1), 112-137. https://doi.org/10.17713/ajs.v54i1.1975
#'
#' @param X A vector representing observed values of the time series.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param meanX The average value of \code{X}. Defaults to \code{mean(X)}.
#' @param x A vector of lag indices. Defaults to the sequence \code{0:length(X)}. Must be at least as large as \code{maxLag + 1}.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the autocovariance/autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the \code{type} parameter.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'standard_est'.}
#' }
#' @export
#'
#' @importFrom stats acf
#'
#' @examples
#' X <- c(1, 2, 3)
#' standard_est(X, pd = FALSE, maxLag = 2, meanX = mean(X))
standard_est <- function(X, pd = TRUE, maxLag = length(X) - 1, x = 0:length(X), type = c("autocovariance", "autocorrelation"), meanX = mean(X)) {
  type <- match.arg(type)
  stopifnot(length(X) > 0, is.numeric(X), is.logical(pd), maxLag >= 0, maxLag <= (length(X) - 1),
            maxLag %% 1 == 0, length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            is.numeric(x), all(!is.na(x)), length(x) > 0)
  retVec <- as.vector(stats::acf(X - meanX, lag.max = maxLag, type = "covariance", plot = FALSE, demean = FALSE)$acf)
  if(!pd) {
    retVec <- retVec * (length(X) / (length(X) - 0:maxLag))
  }
  if(type == 'autocorrelation') {
    retVec <- retVec / retVec[1]
  }

  res <- list(acf=retVec, lags=x[1:(maxLag + 1)], est_type=type, est_used='standard_est')
  return(structure(res, class = "CovEsts"))
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
#' @param estCov A vector whose values are an estimate autocovariance function or a \code{CovEsts} S3 object.
#'
#' @return A numeric vector whose values are an estimate of the semivariogram or an \code{VarioEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{vario}}{A numeric vector whose values are an estimate of the semivariogram.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on, inherited from the argument \code{estCov}.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'to_vario'.}
#' }
#' If a numeric vector is given for the argument \code{estCov}, then a numeric vector output is given, and if a \code{CovEsts} S3 object is given, a \code{VarioEsts} object is given as output.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' estCov <- standard_est(X, meanX=mean(X), maxLag = 2, pd=FALSE)
#' to_vario(estCov)
to_vario <- function(estCov) {
  UseMethod("to_vario")
}

#' @describeIn to_vario Method for numeric vectors.
#' @export
to_vario.default <- function(estCov) {
  stopifnot(length(estCov) > 0, is.numeric(estCov), !any(is.na(estCov)))
  return(estCov[1] - estCov)
}

#' @describeIn to_vario Method for `CovEsts` objects.
#' @export
to_vario.CovEsts <- function(estCov) {
  ret_vec <- to_vario.default(estCov$acf)

  res <- list(vario=ret_vec, lags=estCov$lags, est_used='to_vario')
  return(structure(res, class = "VarioEsts"))
}

#' Computes the Standard Estimator of the Autocovariance Function.
#'
#' This function computes the partial autocorrelation function from an estimated autocovariance or autocorrelation function.
#'
#' @details
#' This function is a translation of the 'uni_pacf' function in src/library/stats/src/pacf.c of the R source code which is an implementation of the Durbin-Levinson algorithm.
#'
#' @param estCov A numeric vector representing an estimated autocovariance or autocorrelation function.
#'
#' @return A vector whose values are an estimate partial autocorrelation function or a \code{CovEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf}}{A numeric vector containing the partial autocorrelation estimates.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on.}
#'  \item{\code{est_type}}{The type of estimate, namely 'partial'.}
#'  \item{\code{est_used}}{The estimator function used, in this case, 'to_pacf'.}
#' }
#' If a numeric vector is given for the argument \code{estCov}, then a numeric vector output is given, and if a \code{CovEsts} S3 object is given, a \code{CovEsts} object is given as output.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' to_pacf(standard_est(X, pd = FALSE, maxLag = 2, meanX = mean(X)))
to_pacf <- function(estCov) {
  UseMethod("to_pacf")
}

#' @describeIn to_pacf Method for `CovEsts` objects.
#' @export
to_pacf.CovEsts <- function(estCov) {
  pacf_vec <- to_pacf.default(estCov$acf)

  res <- list(acf=pacf_vec, lags=1:length(pacf_vec), est_type='partial', est_used='to_pacf')
  return(structure(res, class = "CovEsts"))
}

#' @describeIn to_pacf Method for numeric vectors.
#' @export
to_pacf.default <- function(estCov) {
  stopifnot(length(estCov) > 0, is.numeric(estCov), !any(is.na(estCov)))

  # Force to autocorrelation
  estCov <- estCov / estCov[1]

  # This function is a translation of the 'uni_pacf' function in src/library/stats/src/pacf.c of the R source code.
  n <- length(estCov)
  if(n == 2) {
    return(c(estCov[2]))
  }

  phi <- c(estCov[2], rep(NA, n - 2))
  temp <- rep(NA, n - 1)
  pacf_vec <- numeric(n - 2)
  pacf_vec[1] <- estCov[2]

  for(i in 2:(n - 1)) {
    temp <- phi[1:(i - 1)]
    numer <- estCov[i + 1]
    denom <- estCov[1]

    for(k in 1:(i - 1)) {
      numer <- numer - (phi[k] * estCov[i - k + 1])
      denom <- denom - (phi[k] * estCov[k + 1])
    }

    pacf_val <- numer / denom
    pacf_vec[i] <- pacf_val
    phi[i] <- pacf_val

    if (i == (n - 1)) {
      break
    }

    # update the phi for next run
    # phi_{n, k} <- phi_{n - 1, k} - (phi_{n, n} * phi_{n - 1, n - k})
    for(k in 1:(i - 1)) {
      phi[k] <- temp[k] - (pacf_val * temp[i - k])
    }
  }

  return(pacf_vec)
}

#' Normalise a CovEsts Object
#'
#' @description
#' This function normalises a CovEsts S3 object, that is, turning an autocovariance function into an autocorrelation function.
#' This procedure is a one-way transformation.
#' This function does not support any other type of argument.
#'
#' @param estCov A CovEst S3 object.
#'
#' @return A CovEsts S3 object with a normalised \code{acf} attribute.
#' @export
#'
#' @examples
#' normalise_acf(standard_est(c(1, 2, 3)))
normalise_acf <- function(estCov) {
  UseMethod("normalise_acf")
}

#' @describeIn normalise_acf Method for other objects.
#' @export
normalise_acf.default <- function(estCov) {
  stop("Input must be of class 'CovEsts'. Objects of class ", class(estCov), " are not supported.")
}

#' @describeIn normalise_acf Method for 'CovEsts' objects.
#' @export
normalise_acf.CovEsts <- function(estCov) {
  estCov$acf <- estCov$acf / estCov$acf[1]

  estCov$est_type <- "autocorrelation"

  return(estCov)
}
