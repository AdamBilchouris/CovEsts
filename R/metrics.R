#' Area Between Estimated Autocovariance Functions.
#'
#' This function estimates the area between two estimated autocovariance functions.
#'
#' @details
#' This function estimates the area between two estimated autocovariance functions over a set of lags, from 0 up to \eqn{\tau_{0}.}
#' \deqn{
#' \int_{0}^{\tau_{0}} \left| \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau) \right| d\tau ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' To approximate this integral the trapezoidal rule is used.
#'
#' If \code{lags} is empty, a uniform time grid with a step of 1 will be used which may result in a different area than if lags is specified.
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector representing the second estimated autocovariance function.
#' @param lags An optional vector of lags starting from 0 up until some other lag. If empty, a vector of lags is created starting from 0 until \code{len(est1) - 1}, by 1.
#' @param plot A boolean determining whether a plot should be created. By default, no plot is created.
#'
#' @return A numeric value representing the estimated area between two estimated autocovariance functions.
#' @export
#'
#' @importFrom graphics legend
#' @importFrom graphics lines
#' @importFrom graphics polygon
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2) + rnorm(length(x), sd=0.1)
#' estCov2 <- exp(-x^2.1) + rnorm(length(x), sd=0.1)
#' area_between(estCov1, estCov2, lags=x)
#' area_between(estCov1, estCov2, lags=x, plot = TRUE)
area_between <- function(est1, est2, lags=c(), plot = FALSE) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2), is.logical(plot))

  if(length(lags) != 0) {
    stopifnot(is.vector(lags), length(est1) == length(lags), length(est2) == length(lags))
  }
  # Generate a vector of lags, {0, ..., length(est1) - 1}
  else {
    lags <- seq(0, length(est1) - 1)
  }

  yDiff <- abs(est1 - est2)
  integral_val <- 0
  for(i in 2:length(yDiff)) {
    integral_val_temp <-  ((yDiff[i - 1] + yDiff[i]) / 2) * (lags[i] - lags[i - 1])
    integral_val <- integral_val + integral_val_temp
  }

  if(plot) {
      y_min <- ifelse(min(est1) < min(est2), min(est1), min(est2))
      y_max <- ifelse(max(est1) > max(est2), max(est1), max(est2))
      plot(lags, est1, lwd=2, lty=1, col=2, type='l', ylim=c(y_min, y_max), xlab='Lag (h)', ylab='Estimated Autocovariance')
      graphics::lines(lags, est2, lwd=2, lty=1, col=3)
      # Create a region between the two curves.
      graphics::polygon(c(lags, rev(lags)), c(est2, rev(est1)), col="#AAAAAAAA", border=2, lty=0)
      graphics::legend("topright", inset=0, title="Estimated Autocovariance", c("1", "2", "Area Between"),
             col=c(2, 3, "#AAAAAAAA"), lwd=c(2, 2, NA), lty=c(1, 1, NA), density=c(0, 0, NA), fill=c(NA, NA, "#AAAAAAAA"), border=c(NA, NA, "#000000FF"))
  }

  # \sigma_{i=1}^{n} f(x_{i}) \delta x
  return(integral_val)
}

#' The Maximum Vertical Distance Between Estimated Functions.
#'
#' This function computes the maximum vertical distance between functions.
#'
#' @details
#' This function computes the maximum vertical distance between functions. The vector of function values must be of the same length.
#' \deqn{D(\hat{C}_{1}(\tau), \hat{C}_{2}(\tau)) = \max_{\tau} \left| \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau) \right| ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#' It assumes that the estimated functions are estimated over the same set of lags.
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector representing the second estimated autocovariance function.
#' @param lags An optional vector of lags starting from 0 up until some other lag. If empty, a vector of lags is created starting from 0 until \code{len(est1) - 1}, by 1.
#' @param plot A boolean as to whether a plot should be created. By default, no plot is created.
#'
#' @return A numeric value representing the maximum vertical distance between the two estimated functions.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2) + rnorm(length(x), sd=0.1)
#' estCov2 <- exp(-x^2.1) + rnorm(length(x), sd=0.1)
#' max_distance(estCov1, estCov2)
#' max_distance(estCov1, estCov2, plot = TRUE)
max_distance <- function(est1, est2, lags=c(), plot = FALSE) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2), is.logical(plot))

  if(length(lags) != 0) {
    stopifnot(is.vector(lags), length(est1) == length(lags), length(est2) == length(lags))
  }
  # Generate a vector of lags, {0, ..., length(est1) - 1}
  else {
    lags <- seq(0, length(est1) - 1)
  }

  maxDist <- -Inf

  distVec <- c()

  for(i in 1:length(est1)) {
    d <- abs(est1[i] - est2[i])
    if(is.numeric(d)) {
      distVec <- c(distVec, d)
      if(d > maxDist) {
        maxDist <- d
      }
    }
    else {
      distVec <- c(distVec, 0)
    }
  }

  if(plot) {
    plot(lags, distVec, type='o', xlab='Lag (h)', ylab='Distance')
  }

  return(maxDist)
}

#' Create a Cyclic Matrix for a Given Vector.
#'
#' This helper function creates a cyclic matrix from a given vector \eqn{v}.
#'
#' @details
#' This function creates a cylic matrix for a given vector \eqn{v}.
#' If \eqn{v = \{v_{0}, v_{1} , \dots , v_{N-1} , v_{N} \},} then the cyclic matrix will has the form
#' \deqn{\left[ {\begin{array}{ccccc}
#' v_{0}     & v_{1}     & \cdots & v_{N - 1} & v_{N}     \\
#' v_{1}     & v_{0}     & \cdots & v_{N - 2} & v_{N - 1} \\
#' \vdots    & \vdots    & \ddots & \vdots    & \vdots    \\
#' v_{N - 1} & v_{N- 2}  & \cdots & v_{0}     & v_{1}     \\
#' v_{N}     & v_{N - 1} & \cdots & v_{1}     & v_{0}     \\
#' \end{array}} \right]
#' }
#'
#' @param v A numeric vector.
#'
#' @return A cyclic matrix.
#' @export
#'
#' @examples
#' v <- c(1, 2, 3)
#' create_cyclic_matrix(v)
create_cyclic_matrix <- function(v) {
 stopifnot(is.numeric(v), !any(is.na(v)), length(v) >= 1)
  n <- length(v)
  mat <- matrix(NA, nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      mat[i, j] <- v[1 + abs(i - j)]
    }
  }
  return(mat)
}

#' Compute the Spectral Norm Between Estimated Functions.
#'
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#' This function is intended for estimates over lags with a constant difference.
#'
#' @details
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#' Let \eqn{D(\tau) = \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau),}
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' A matrix \eqn{D} is created from \eqn{D(\cdot)},
#' \deqn{\left[ {\begin{array}{ccccc}
#' D(0)            & D(\tau_{1})     & \cdots & D(\tau_{N - 1}) & D(\tau_{N})     \\
#' D(\tau_{1})     & D(0)            & \cdots & D(\tau_{N - 2}) & D(\tau_{N - 1}) \\
#' \vdots          & \vdots          & \ddots & \vdots          & \vdots       \\
#' D(\tau_{N - 1}) & D(\tau_{N - 2}) & \cdots & D(0)            & D(\tau_{1})     \\
#' D(\tau_{N})     & D(\tau_{N - 1}) & \cdots & D(\tau_{1})     & D(0)         \\
#' \end{array}} \right] ,
#' }
#' over a set of lags \eqn{\{0, \tau_{1}, \dots , \tau_{N} \}.}
#' This matrix is created in [create_cyclic_matrix].
#'
#' The spectral norm is simply the largest eigenvalue of \eqn{D.}
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector representing the second estimated autocovariance function.
#'
#' @return The spectral norm of the differences between the two functions.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2) + rnorm(length(x), sd=0.1)
#' estCov2 <- exp(-x^2.1) + rnorm(length(x), sd=0.1)
#' spectral_norm(estCov1, estCov2)
spectral_norm <- function(est1, est2) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2))
  det1 <- norm(create_cyclic_matrix((est1 - est2)), type = "2")
  if(is.numeric(det1)) {
    return(det1)
  }
  return(NA)
}

#' Check if an Autocovariance Function Estimate is Positive-Definite or Not.
#'
#' This function checks if an autocovariance function estimate is positive-definite or not by determining if the eigenvalues of the corresponding matrix (see the Details section) are all positive.
#'
#' @details
#' For an autocovariance function estimate \eqn{\hat{C}(\cdot)} over a set of lags separated by a constant difference \eqn{\{0, h_{1} , h_{2} , \dots , h_{n} \},}
#' construct the symmetric matrix
#' \deqn{\left[ {\begin{array}{ccccc}
#' \hat{C}(0)         & \hat{C}(h_{1})     & \cdots & \hat{C}(h_{N - 1}) & \hat{C}(h_{N})     \\
#' \hat{C}(h_{1})     & \hat{C}(0)         & \cdots & \hat{C}(h_{N - 2}) & \hat{C}(h_{N - 1}) \\
#' \vdots             & \vdots             & \ddots & \vdots             & \vdots             \\
#' \hat{C}(h_{N - 1}) & \hat{C}(h_{N - 2}) & \cdots & \hat{C}(0)         & \hat{C}(h_{1})     \\
#' \hat{C}(h_{N})     & \hat{C}(h_{N - 1}) & \cdots & \hat{C}(h_{1})     & \hat{C}(0)         \\
#' \end{array}} \right] .
#' }
#'
#' The eigendecomposition of this matrix is computed to determine if all eigenvalues are positive, if so, the estimated autocovariance function is assumed to be positive-definite.
#'
#' @param est A numeric vector representing an estimated autocovariance function.
#'
#' @return A boolean where \code{TRUE} denotes a positive-definite autocovariance function estimate and \code{FALSE} for an estimate that is not positive-definite.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov <- exp(-x^2)
#' check_pd(estCov)
check_pd <- function(est) {
  stopifnot(is.numeric(est), length(est) >= 1, !any(is.na(est)))
  cyclic_mat <- create_cyclic_matrix(est)
  eigens <- zapsmall(eigen(cyclic_mat)$values)
  return(!any(eigens < 0) && !all(eigens == 0))
}
