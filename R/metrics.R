#' Area Between Estimated Autocovariance Functions.
#'
#' This function estimates the area between two estimated autocovariance functions.
#'
#' @details
#' This function estimates the area between two estimated autocovariance functions over a set of lags, from 0 up to \eqn{h_{n}} defined by
#' \deqn{
#' \int_{0}^{h_{n}} \left| \hat{C}_{1}(h) - \hat{C}_{2}(h) \right| dh ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' To approximate this integral the trapezoidal rule is used.
#'
#' If \code{lags} is empty, a uniform time grid with a step of 1 will be used which may result in a different area than if lags is specified.
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector of the same length as \code{est1} representing the second estimated autocovariance function
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
#' estCov1 <- exp(-x^2)
#' estCov2 <- exp(-x^2.1)
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

  return(integral_val)
}

#' Maximum Vertical Distance Between Estimated Functions.
#'
#' This function computes the maximum vertical distance between functions.
#'
#' @details
#' This function computes the maximum vertical distance between functions:
#' \deqn{D(\hat{C}_{1}(h), \hat{C}_{2}(h)) = \displaystyle \max_{h} \left| \hat{C}_{1}(h) - \hat{C}_{2}(h) \right| ,
#' }
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#' It assumes that the estimated values are given for the same set of lags.
#' The vectors of the function values must be of the same length.
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector of the same length as \code{est1} representing the second estimated autocovariance function
#' @param lags An optional vector of lags starting from 0 up until some other lag. If empty, a vector of lags is created starting from 0 until \code{len(est1) - 1}, by 1.
#' @param plot A boolean as to whether a plot should be created. By default, no plot is created.
#'
#' @return A numeric value representing the maximum vertical distance between the two estimated functions.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2)
#' estCov2 <- exp(-x^2.1)
#' max_distance(estCov1, estCov2, lags=x)
#' max_distance(estCov1, estCov2, lags=x, plot = TRUE)
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
    plot(lags, distVec, type='l', xlab='Lag (h)', ylab='Distance')
  }

  return(maxDist)
}

#' Create a Cyclic Matrix for a Given Vector.
#'
#' This helper function creates a symmetric matrix from a given vector \eqn{v}.
#'
#' @details
#' This function creates a symmetric matrix for a given vector \eqn{v}.
#' If \eqn{v = \{v_{0}, v_{1} , \dots , v_{N-1} , v_{N} \},} then the symmetric matrix will has the form
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
#' @return A symmetric matrix.
#' @export
#'
#' @importFrom stats toeplitz
#'
#' @examples
#' v <- c(1, 2, 3)
#' cyclic_matrix(v)
cyclic_matrix <- function(v) {
  stopifnot(is.numeric(v), !any(is.na(v)), length(v) >= 1)
  return(toeplitz(v))
}

#' Compute the Spectral Norm Between Estimated Functions.
#'
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#' This function is intended for estimates over lags with a constant difference.
#'
#' @details
#' This function computes the spectral norm of the difference of two estimated autocovariance functions.
#' Let \eqn{D(h) = \hat{C}_{1}(h) - \hat{C}_{2}(h),}
#' where \eqn{\hat{C}_{1}(\cdot)} and \eqn{\hat{C}_{2}(\cdot)} are estimated autocovariance functions.
#'
#' A matrix \eqn{D} is created from \eqn{D(\cdot)},
#' \deqn{\left[ {\begin{array}{ccccc}
#' D(h_{0})     & D(h_{1})     & \cdots & D(h_{n - 1}) & D(h_{n})     \\
#' D(h_{1})     & D(h_{0})     & \cdots & D(h_{n - 2}) & D(h_{n - 1}) \\
#' \vdots       & \vdots       & \ddots & \vdots       & \vdots       \\
#' D(h_{n - 1}) & D(h_{n - 2}) & \cdots & D(h_{0})     & D(h_{1})     \\
#' D(h_{n})     & D(h_{n - 1}) & \cdots & D(h_{1})     & D(h_{0})     \\
#' \end{array}} \right] ,
#' }
#' over a set of lags \eqn{\{h_{0}, h_{1}, \dots , h_{N} \}.}
#' This matrix is created by [cyclic_matrix].
#'
#' The spectral norm is defined as the largest eigenvalue of \eqn{D.}
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector of the same length as \code{est1} representing the second estimated autocovariance function
#'
#' @return The spectral norm of the differences between the two functions.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2)
#' estCov2 <- exp(-x^2.1)
#' spectral_norm(estCov1, estCov2)
spectral_norm <- function(est1, est2) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2))
  det1 <- norm(cyclic_matrix((est1 - est2)), type = "2")
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
#' For an autocovariance function estimate \eqn{\hat{C}(\cdot)} over a set of lags separated by a constant difference \eqn{\{h_{0}, h_{1} , h_{2} , \dots , h_{n} \},}
#' construct the symmetric matrix
#' \deqn{\left[ {\begin{array}{ccccc}
#' \hat{C}(h_{0})     & \hat{C}(h_{1})     & \cdots & \hat{C}(h_{n - 1}) & \hat{C}(h_{n})     \\
#' \hat{C}(h_{1})     & \hat{C}(h_{0})     & \cdots & \hat{C}(h_{n - 2}) & \hat{C}(h_{n - 1}) \\
#' \vdots             & \vdots             & \ddots & \vdots             & \vdots             \\
#' \hat{C}(h_{n - 1}) & \hat{C}(h_{n - 2}) & \cdots & \hat{C}(h_{0})     & \hat{C}(h_{1})     \\
#' \hat{C}(h_{n})     & \hat{C}(h_{n - 1}) & \cdots & \hat{C}(h_{1})     & \hat{C}(h_{0})     \\
#' \end{array}} \right] .
#' }
#'
#' The eigendecomposition of this matrix is computed to determine if all eigenvalues are positive. If so, the estimated autocovariance function is assumed to be positive-definite.
#'
#' @param est A numeric vector or corresponding cyclic matrix representing an estimated autocovariance function.
#'
#' @return A boolean where \code{TRUE} denotes a positive-definite autocovariance function estimate and \code{FALSE} for an estimate that is not positive-definite.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov <- exp(-x^2)
#' check_pd(estCov)
#' check_pd(cyclic_matrix(estCov))
check_pd <- function(est) {
  stopifnot(is.numeric(est), length(est) >= 1, !any(is.na(est)))

  if(is.matrix(est)) {
    eigens <- zapsmall(eigen(est)$values)
    return(!any(eigens < 0) && !all(eigens == 0))
  }

  cyclic_mat <- cyclic_matrix(est)
  eigens <- zapsmall(eigen(cyclic_mat)$values)
  return(!any(eigens < 0) && !all(eigens == 0))
}

#' MSE Between Estimated Autocovariance Functions.
#'
#' This function computes the mean-square difference/error between two autocovariance functions (estimated or theoretical).
#'
#' @details
#' This function computes the mean-square difference/error (MSE) between two estimated autocovariance functions (estimated or theoretical).
#' The MSE is defined as
#' \deqn{\frac{1}{n} \sum_{i=0}^{n} \left(\widehat{C}_{1}(h_{i}) - \widehat{C}_{2}(h_{i})\right)^{2}}
#' over a set of lags \eqn{\{h_{0}, h_{1} , h_{2} , \dots , h_{n} \}.}
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector of the same length as \code{est1} representing the second estimated (or theoretical) autocovariance function
#'
#' @return A numeric value representing the MSE between two autocovariance functions (estimated or theoretical).
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2)
#' estCov2 <- exp(-x^2.1)
#' mse(estCov1, estCov2)
mse <- function(est1, est2) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2))
  return(mean((est1 - est2)^2))
}

#' Hilbert-Schmidt Norm Between Estimated Autocovariance Functions.
#'
#' This function computes the Hilbert-Schidmt norm between two estimated autocovariance functions.
#'
#' @details
#' This function computes the Hilbert-Schidmt norm between two estimated autocovariance functions.
#' The Hilbert-Schmidt norm of a matrix
#' \deqn{D = \left[(d_{i,j})_{1 \le i,j \le n}\right] = \left[ {\begin{array}{ccccc}
#' D(h_{0})     & D(h_{1})     & \cdots & D(h_{n - 1}) & D(h_{n})     \\
#' D(h_{1})     & D(h_{0})     & \cdots & D(h_{n - 2}) & D(h_{n - 1}) \\
#' \vdots       & \vdots       & \ddots & \vdots       & \vdots       \\
#' D(h_{n - 1}) & D(h_{n - 2}) & \cdots & D(h_{0})     & D(h_{1})     \\
#' D(h_{n})     & D(h_{n - 1}) & \cdots & D(h_{1})     & D(h_{0})     \\
#' \end{array}} \right] ,
#' }
#' over a set of lags \eqn{\{h_{0}, h_{1}, \dots , h_{N} \},} where \eqn{D(h) = \hat{C}_{1}(h) - \hat{C}_{2}(h),}
#' is defined as \deqn{{\left\Vert D \right\Vert}_{HS} = \sqrt{\sum_{i,j} d_{i, j}^{2}}.}
#'
#' @param est1 A numeric vector representing the first estimated autocovariance function.
#' @param est2 A numeric vector of the same length as \code{est1} representing the second estimated autocovariance function
#'
#' @return A numeric value representing the estimated Hilbert-Schmidt norm between two estimated autocovariance functions.
#' @export
#'
#' @examples
#' x <- seq(0, 5, by=0.1)
#' estCov1 <- exp(-x^2)
#' estCov2 <- exp(-x^2.1)
#' hilbert_schmidt(estCov1, estCov2)
hilbert_schmidt <- function(est1, est2) {
  stopifnot(is.numeric(est1), is.numeric(est2), length(est1) >= 1, length(est2) >= 1, !any(is.na(est1)), !any(is.na(est2)),
            length(est1) == length(est2))
  return(sqrt(sum(cyclic_matrix(est1 - est2)^2)))
}

#' Compute the Nearest Positive-Definite Matrix.
#'
#' This function computes the nearest positive-definite matrix to some matrix \eqn{A}.
#'
#' @details
#' This function computes the nearest positive-definite matrix to some matrix \eqn{A}.
#'
#' The procedure to do so is as follows
#' For a matrix \eqn{X}, compute the symmetric matrix \eqn{B = (A + A^{T}) / 2.}
#' Let \eqn{B = UH} be the polar decomposition of B.
#' The nearest positive-definite matrix to \eqn{X} is \eqn{X_{F} = (B + H) / 2.}
#'
#' Unlike [shrinking], only an autocorrelation matrix can be returned, not an autocovariance function.
#'
#' The implementation is a translation of https://au.mathworks.com/matlabcentral/fileexchange/42885-nearestspd#functions_tab .
#'
#' @references
#' Higham, N. J. (1988). Computing a nearest symmetric positive semidefinite matrix. Linear Algebra and its Applications, 103, 103–118. 10.1016/0024-3795(88)90223-6
#'
#' D'Errico, J. (2025). nearestSPD (https://www.mathwor ks.com/matlabcentral/fileexchange/42885-nearestspd), MATLAB Central File Exchange. Retrieved August 2, 2025.
#'
#' @param X Either a numeric vector or a square matrix. If a vector is provided, a matrix will be created of the form found in [cyclic_matrix].
#'
#' @return The closest positive-definite autocorrelation matrix.
#' @export
#'
#' @examples
#' X <- c(1, 0, -1.1)
#' nearest_pd(X)
#' check_pd(nearest_pd(X))
nearest_pd <- function(X) {
  stopifnot(is.numeric(X), !any(is.na(X)), length(X) > 0)
  if(!is.matrix(X)) {
    X <- cyclic_matrix(X)
  }
  # https://au.mathworks.com/matlabcentral/fileexchange/42885-nearestspd#functions_tab
  B <- (X + t(X)) / 2

  svd_mat <- svd(X)
  # https://en.wikipedia.org/wiki/Polar_decomposition#Relation_to_the_SVD
  H <- svd_mat$v %*% diag(svd_mat$d) %*% t(svd_mat$v)

  X_F <- (B + H) / 2
  X_F <- (X_F + t(X_F)) / 2

  ident_mat <- diag(nrow(X_F))
  k <- 1
  while(!all(zapsmall(eigen(X_F)$values) >= 0)) {
    min_eig <- min(eigen(X_F)$values)
    X_F <- X_F + (ident_mat * (min_eig * k^2 + abs(min_eig - (min_eig + .Machine$double.eps * min_eig))))
  }
  return(X_F)
}
