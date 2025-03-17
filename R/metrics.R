#' Area between estimated functions.
#'
#' This function computes the area between functions.
#'
#' @details
#' This function computes the area between the two estimated functions over a set of lags, from 0 up to \eqn{\tau_{0}.}
#' \deqn{
#' \int_{0}^{\tau_{0}} \left| \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau) \right| d\tau ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' If \code{lags} is empty a difference of 1 will be used which results in a different area than if lags is specified.
#'
#' @param est1 A numeric vector representing the first estimated covariance function.
#' @param est2 A numeric vector representing the second estimated covariance function.
#' @param lags An optional vector of lags starting from 0 up until some other lag. If empty, a vector of lags is created starting from 0 until \code{len(est1) - 1}, by 1.
#' @param plot A boolean as to whether a plot should be created. By default, no plot is created.
#'
#' @return A numeric value representing the area between the estimated functions.
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
  # (b-a) / n
  weight <- (lags[length(lags)] - lags[1]) / length(lags)

  if(plot) {
      y_min <- ifelse(min(est1) < min(est2), min(est1), min(est2))
      y_max <- ifelse(max(est1) > max(est2), max(est1), max(est2))
      plot(lags, est1, lwd=2, lty=1, col=2, type='l', ylim=c(y_min, y_max))
      graphics::lines(lags, est2, lwd=2, lty=1, col=3)
      # Create a region between the two curves.
      graphics::polygon(c(lags, rev(lags)), c(est2, rev(est1)), col="#AAAAAAAA", border=2, lty=0)
      graphics::legend("topright", inset=0, title="Estimated Autocovariance", c("1", "2", "Area Between"),
             col=c(2, 3, "#AAAAAAAA"), lwd=c(2, 2, NA), lty=c(1, 1, NA), density=c(0, 0, NA), fill=c(NA, NA, "#AAAAAAAA"), border=c(NA, NA, "#000000FF"))
  }

  # \sigma_{i=1}^{n} f(x_{i}) \delta x
  return(sum(yDiff * weight))
}

#' The maximum vertical distance between estimated functions.
#'
#' This function computes the maximum vertical distance between functions.
#'
#' @details
#' This function computes the maximum vertical distance between functions. The two functions must be of the same length.
#' \deqn{D(\hat{C}_{1}(\tau), \hat{C}_{2}(\tau)) = \max_{\tau} \left| \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau) \right| ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#' It assumes that the estimated functions are estimated over the same set of lags.
#'
#' @param est1 A numeric vector representing the first estimated covariance function.
#' @param est2 A numeric vector representing the second estimated covariance function.
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
max_distance <- function(est1, est2, plot = FALSE) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2), is.logical(plot))
  maxDist <- -Inf

  # Compute the distances, I doubt it takes up much memory, it is being computed anyway
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
    plot(distVec, type='o')
  }

  return(maxDist)
}

#' Create a cyclic matrix for a given vector.
#'
#' This helper function creates a cyclic matrix for a given vector \eqn{v}.
#'
#' @details
#' This function creates a cylic matrix for a given vector \eqn{v}.
#' If \eqn{v = \{v_{0}, v_{1} , \dots , v_{N-1} , v_{N} \},} then the cyclic matrix will
#' be of the form
#' \deqn{\left[ {\begin{array}{ccccc}
#' v_{0}     & v_{1}     & \cdots & v_{N - 1} & v_{N}     \\
#' v_{1}     & v_{0}     & \cdots & v_{N - 2} & v_{N - 1} \\
#' \vdots    & \vdots    & \ddots & \vdots    & \vdots    \\
#' v_{N - 1} & v_{N- 2}  & \cdots & v_{0}     & v_{1}     \\
#' v_{N}     & v_{N - 1} & \cdots & v_{1}     & v_{0}     \\
#' \end{array}} \right] ,
#' }
#'
#' @param v A vector of values.
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

#' Compute spectral norm between estimated functions.
#'
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#'
#' @details
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#' Let \eqn{\hat{D}(\tau) = \hat{C}_{1}(\tau) - \hat{C}_{2}(\tau),}
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' A matrix \eqn{D} is created from \eqn{\hat{D}(\cdot)},
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
#' The spectral norm is simply the largest singular value of \eqn{D.}
#'
#' @param est1 A numeric vector representing the first estimated covariance function.
#' @param est2 A numeric vector representing the second estimated covariance function.
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
