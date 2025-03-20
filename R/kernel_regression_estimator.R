#' Compute \eqn{X_{ij}} Matrix
#'
#' This helper function computes the matrix of pairwise values, \eqn{X_{ij},} for the kernel regression estimator,
#' \deqn{X_{ij} = (X_{i} - \bar{X}) (X_{j} - \bar{X}) .}
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#'
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the nonparametric estimation of covariance functions. The Annals of Statistics (Vol. 22, Issue 4, pp. 2115–2134). 10.1214/aos/1176325774
#'
#' @param X A vector of values.
#'
#' @return A matrix of size \eqn{N \times N}, where \eqn{N} is the length of the vector X.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' Xij_mat(X)
Xij_mat <- function(X) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)))

  meanX <- mean(X)
  Xij <- matrix(nrow=length(X), ncol=length(X))
  for(i in 1:length(X)) {
    for(j in 1:length(X)) {
      Xij[i, j] <- (X[i] - meanX) * (X[j] - meanX)
    }
  }

  return(Xij)
}

#' Compute \eqn{\rho(T_{1})} for the truncated kernel regression estimator.
#'
#' This helper function computes \eqn{\rho(T_{1})} for the truncated kernel regression estimator.
#'
#' @details
#' This function computes the following value,
#' \deqn{
#' \hat{\rho}(T_{1}) = \left( \sum_{i} \sum_{j} \check{X}_{ij} K((T_{1} - t_{ij}) / h) \right) \left( \sum_{i} \sum_{j} K((T_{1} - t_{ij}) / h) \right)^{-1},
#' }
#' where \eqn{\check{X}_{ij} = (X(t_{i}) - \bar{X}) (X(t_{j}) - \bar{X})}, and \eqn{t_{ij} = t_{i} - t_{j}.}
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of Nonparametric Estimators of Autocovariance for Stationary Random Fields. Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#'
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the nonparametric estimation of covariance functions. The Annals of Statistics (Vol. 22, Issue 4, pp. 2115–2134). 10.1214/aos/1176325774
#'
#' @param x A vector of lags.
#' @param meanX The average value of X.
#' @param T1 The first trunctation point.
#' @param h Bandwidth parameter.
#' @param xij_mat The matrix of pairwise covariance values.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j". Alternatively, a custom kernel function can be provided, see [compute_corrected_standard_est]'s example.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel_symm] for parameters.
#' @param custom_kernel If a custom kernel is to be used or not.
#'
#' @return The estimated covariance function at \eqn{T_{1}}.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' rho_T1(1:4, mean(X), 1, 0.1, Xij_mat(X), "gaussian", c(), FALSE)
rho_T1 <- function(x, meanX, T1, h, xij_mat, kernel_name="gaussian", kernel_params=c(), custom_kernel = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            length(T1) == 1, is.numeric(T1), !is.na(T1), T1 > 0, length(h) == 1, is.numeric(h), !is.na(h), h > 0,
            is.numeric(xij_mat), is.matrix(xij_mat), !any(is.na(xij_mat)), is.logical(custom_kernel))

  numerators <- c()
  denominators <- c()

  # Duplicating code as to only have 1 if statement, not length(x) if statements.
  if(custom_kernel) {
    for(i in 1:length(x)) {
      tij <- x[i] - x
      t_tij <- T1 - tij
      X_ij <- xij_mat[i, ]
      K_ij <- get(kernel_name)(t_tij, h, kernel_params)
      numerator <- K_ij * X_ij
      denominator <- K_ij

      numerators <- c(numerators, sum(numerator))
      denominators <- c(denominators, sum(denominator))
    }
  }
  else {
    stopifnot(kernel_name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    for(i in 1:length(x)) {
      tij <- x[i] - x
      t_tij <- T1 - tij
      X_ij <- xij_mat[i, ]
      K_ij <- get("kernel_symm")(t_tij, kernel_name, c(h, kernel_params[1], kernel_params[2]))
      numerator <- K_ij * X_ij
      denominator <- K_ij

      numerators <- c(numerators, sum(numerator))
      denominators <- c(denominators, sum(denominator))
    }
  }

  return( ( sum(numerators) / sum(denominators) ) )
}

#' Compute 1D Discrete Cosine Transform
#'
#' This helper function computes the Type-II discrete cosine transform.
#'
#' The Type-II discrete cosine transform is obtained using [stats::fft]. Note this method does not utilise Makhoul's method for computing the DCT.
#' A new input \eqn{y} is created from \eqn{X} where \eqn{y_{2n} = 0}, \eqn{y_{2n+1} = X_{n}} for \eqn{0 \le n < N}, \eqn{y_{2N} = 0} and \eqn{y_{4N - n} = y_{n}}
#' for \eqn{0 < n < 2N}, where N is the length of \eqn{X}.
#' After this, the is done (in code): \code{dct <- 0.5 * Re(stats::fft(Y))[1:(length(Y) / 4)]}, which gives the Type-II discrete cosine transform.
#'
#' @references
#' Ochoa-Dominguez, H., & Rao, K.R. (2019). Discrete Cosine Transform, Second Edition. CRC Press. 10.1201/9780203729854
#'
#' @param X A vector of values for which the discrete cosine transform is being computed.
#'
#' @return A vector of discrete cosine transform values.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' dct_1d(X)
dct_1d <- function(X) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)))
  newX <- c(X, rev(X))
  zerosX <- rep(0, length(newX))
  newX <- c(rbind(zerosX, newX))

  dct_seq <- 0.5*Re(stats::fft(newX))[1:(length(newX)/4)]
  return(dct_seq)
}

#' Compute 1D Inverse Discrete Cosine Transform
#'
#' This helper function computes the inverse of the Type-II discrete cosine transform.
#'
#' @details
#' The Type-II inverse discrete cosine transform is computed using [stats::fft].
#' The original spectrum is obtained from the input X and then an inverse FFT is applied to obtain
#' the estimate.
#'
#' The original spectrum, \code{dct_full}, from \code{X} is obtained as follows:
#' \code{dct_full <- c(X, 0, -X[-1], 0, rev(X[-1]))}. After this, an inverse FFT is applied,
#' \code{idct <- Re(stats::fft(dct_full, inverse = T)) * (2 / length(dct_full))}, and then the first half of inverse transform is taken, and the zeroes are dropped,
#' which gives the untransformed X.
#'
#' @references
#' Ochoa-Dominguez, H., & Rao, K.R. (2019). Discrete Cosine Transform, Second Edition. CRC Press. 10.1201/9780203729854
#'
#' endolith (2013, September 4). Fast Cosine Transform via FFT. https://dsp.stackexchange.com/a/10606 (archive: https://web.archive.org/web/20240423074416/https://dsp.stackexchange.com/questions/2807/fast-cosine-transform-via-fft).
#'
#' @param X A vector of values for which the discrete cosine transform is being computed.
#'
#' @return A vector with the inverse transform values.
#' @export
#'
#' @examples
#' X <- dct_1d(c(1, 2, 3))
#' idct_1d(X)
idct_1d <- function(X) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)))
  # Reconstruct original frequency (i.e. the full dct_seq from dct_1d)
  dct_full <- c(X, 0, -rev(X), -X[-1], 0, rev(X[-1]))

  # Perform the inversion
  idct_vals <- Re(stats::fft(dct_full, inverse = T)) * (2 / length(dct_full))

  # Recover original signal, drop 0s
  seqEven <- seq(2, 2*length(X), by=2)
  idct_vals_2 <- idct_vals[seqEven]

  return(idct_vals_2)
}

#' Compute the Truncated Kernel Regression Estimator.
#'
#' This function computes the truncated kernel regression estimator,
#' \deqn{\hat{\rho}_{1}(t)  = \left\{ \begin{array}{ll}
#' \hat{\rho}(t) & 0 \leq t \leq T_{1} \\
#' \hat{\rho}(T_{1}) (T_{2} - t)(T_{2} - T_{1})^{-1} & T_{1} < t \leq T_{2} \\
#' 0 & t > T_{2}
#' \end{array} , \right. }
#' where \eqn{\hat{\rho}(\cdot)} is the kernel regression estimator, see [compute_adjusted_est].
#'
#' @details
#' This function computes the truncated kernel regression estimator,
#' \deqn{\hat{\rho}_{1}(t) = \left\{ \begin{array}{ll}
#' \hat{\rho}(t) & 0 \leq t \leq T_{1} \\
#' \hat{\rho}(T_{1}) (T_{2} - t)(T_{2} - T_{1})^{-1} & T_{1} < t \leq T_{2} \\
#' 0 & t > T_{2}
#' \end{array} \right. }
#' where \eqn{\hat{\rho}(\cdot)} is the kernel regression estimator, see [compute_adjusted_est].
#'
#' Compared to [compute_adjusted_est], this function brings down the estimate to zero linearly between \eqn{T_{1}} and \eqn{T_{2}}.
#' In the case of short-range dependence, this may be beneficial as it can remove estimation artefacts at large lags.
#'
#' To make this estimator positive-definite, the following procedure is used:
#' 1. Take the discrete cosine transform
#' \eqn{\mathcal{F}^{c}(\hat{\rho}_{1}(t))}.
#' 2. Find the smallest frequency where its associated value in the spectral domain is negative
#' \deqn{\hat{\theta} = \inf\{ \theta > 0 :  \mathcal{F}^{c}(\hat{\rho}_{1}(t)) < 0\}.}
#' 3. Set all values starting at the frequency to zero.
#' 4. Perform the inversion.
#'
#' This ensures the covariance function estimate is positive-definite. If \eqn{\hat{\theta}} is the first nonzero sample frequency,
#' the entire spectrum, apart from an impulse at zero, is zero, resulting in an adjusted function that is a horizontal line whose value is the area of the of the estimated function prior to the adjustment.
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#'
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the nonparametric estimation of covariance functions. The Annals of Statistics (Vol. 22, Issue 4, pp. 2115–2134). 10.1214/aos/1176325774
#'
#' @param X A vector representing the values of the process.
#' @param x A vector of lags.
#' @param t The values at which the covariance function is calculated at.
#' @param T1 The first truncation point, \eqn{T_{1} > 0.}
#' @param T2 The second truncation point, \eqn{T_{2} > T_{1} > 0.}
#' @param h Bandwidth parameter.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j". Alternatively, a custom kernel function can be provided, see [compute_corrected_standard_est]'s example.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel_symm] for parameters.
#' @param custom_kernel If a custom kernel is to be used or not.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param type Compute either the 'covariance' or 'correlation'. Defaults to 'covariance'.
#' @param meanX The average value of X. Defaults to \code{mean(X)}.
#'
#' @return A vector whose values are the truncated kernel regression estimator.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_truncated_est(X, 1:4, 1:3, 1, 2, 0.1,
#'                   "gaussian", c(), FALSE, TRUE, meanX = mean(X))
compute_truncated_est <- function(X, x, t, T1, T2, h, kernel_name="gaussian",  kernel_params=c(), custom_kernel = FALSE, pd = TRUE, type='covariance', meanX = mean(X)) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), length(x) >= 1, is.numeric(x), !any(is.na(x)),
            length(meanX) == 1, is.numeric(meanX), !is.na(meanX), is.numeric(t), length(t) >= 1,
            length(T1) == 1, is.numeric(T1), !is.na(T1), T1 > 0, length(T2) == 1, is.numeric(T2),
            !is.na(T2), T2 > T1, length(h) == 1, is.numeric(h), h > 0, is.logical(custom_kernel),
            type %in% c('covariance', 'correlation'), is.logical(pd))

  xij_mat <- Xij_mat(X)
  rhoT1 <- rho_T1(x, meanX, T1, h, xij_mat, kernel_name, kernel_params, custom_kernel)

  outer_x_x <- outer(x, x, '-')
  vals_truncated_1 <- c()

  t_vals <- t[t <= T1]
  if(custom_kernel) {
    for(tt in t_vals) {
      K_ij <- get(kernel_name)(tt - outer_x_x, h, kernel_params)
      vals_truncated_1 <- c(vals_truncated_1, sum(K_ij * xij_mat) / sum(K_ij))
    }
  }
  else {
    stopifnot(kernel_name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    for(tt in t_vals) {
      K_ij <- kernel_symm(tt - outer_x_x, kernel_name, c(h, kernel_params[1], kernel_params[2]))
      vals_truncated_1 <- c(vals_truncated_1, sum(K_ij * xij_mat) / sum(K_ij))
    }
  }

  vals_truncated_1 <- c(vals_truncated_1, rhoT1 * (T2 - t[t > T1 & t <= T2]) * (T2 - T1)^(-1))
  vals_truncated_1 <- c(vals_truncated_1, rep(0, length(t[t > T2])))

  if(pd) {
    # Perform DCT
    vals_truncated_1_dct <- dct_1d(vals_truncated_1)

    # Find the first nonzero negative index.
    firstNeg <- which(vals_truncated_1_dct[-1] < 0)[1] + 1

    if(length(firstNeg) == 0 || is.na(firstNeg)) {
      return(vals_truncated_1)
    }

    # Set values after the first negative index onwards
    vals_truncated_1_dct[firstNeg:length(vals_truncated_1_dct)] <- 0

    # Inversion
    vals_truncated_1_idct <- idct_1d(vals_truncated_1_dct)

    if(type == 'correlation') {
      vals_truncated_1_idct <- vals_truncated_1_idct / vals_truncated_1_idct[1]
    }

    return(vals_truncated_1_idct)
  }

  if(type == 'correlation') {
    vals_truncated_1 <- vals_truncated_1 / vals_truncated_1[1]
  }

  return(vals_truncated_1)
}

#' Compute the Kernel Regression Estimator.
#'
#' The kernel regression estimator is defined as
#' \deqn{
#' \hat{\rho}(t) = \left( \sum_{i} \sum_{j} \check{X}_{ij} K((t - t_{ij}) / h) \right) \left( \sum_{i} \sum_{j} K((t - t_{ij}) / h) \right)^{-1},
#' }
#' where \eqn{\check{X}_{ij} = (X(t_{i}) - \bar{X}) (X(t_{j}) - \bar{X})}, and \eqn{t_{ij} = t_{i} - t_{j}.}
#'
#' @details
#' The kernel regression estimator of an autocovariance function is defined as
#' \deqn{
#' \hat{\rho}(t) = \left( \sum_{i} \sum_{j} \check{X}_{ij} K((t - t_{ij}) / h) \right) \left( \sum_{i} \sum_{j} K((t - t_{ij}) / h) \right)^{-1},
#' }
#' where \eqn{\check{X}_{ij} = (X(t_{i}) - \bar{X}) (X(t_{j}) - \bar{X})}, and \eqn{t_{ij} = t_{i} - t_{j}.}
#'
#' If \code{pd} is \code{TRUE} estimator will be made positive-definite through the following procedure
#' 1. Take the discrete cosine transform,
#' \eqn{\mathcal{F}^{c}(\hat{\rho}(t))}.
#' 2. Set all negative values to zero,
#' \eqn{\hat{\mathcal{F}}^{c}(\hat{\rho}(t)) = \mathcal{F}(\hat{\rho}(t))} if \eqn{\mathcal{F}(\hat{\rho}(t)) > 0} and \eqn{\hat{\mathcal{F}}^{c}(\hat{\rho}(t)) = 0} otherwise,
#' for all sample frequencies.
#' 3. Perform the inversion to obtain a new estimator.
#'
#' This ensures the covariance estimate is positive-definite.
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#'
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the nonparametric estimation of covariance functions. The Annals of Statistics (Vol. 22, Issue 4, pp. 2115–2134). 10.1214/aos/1176325774
#'
#' @param X A vector representing the values of the process.
#' @param x A vector of lags.
# #' @param meanX The average value of X.
#' @param t The values at which the covariance function is calculated at.
#' @param h Bandwidth parameter.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j". Alternatively, a custom kernel function can be provided, see [compute_corrected_standard_est]'s example.
#' @param kernel_params A vector of parameters of the kernel function. See [kernel_symm] for parameters.
#' @param custom_kernel If a custom kernel is to be used or not.
#' @param pd Whether a positive-definite estimate should be used. Defaults to \code{TRUE}.
#' @param type Compute either the 'covariance' or 'correlation'. Defaults to 'covariance'.
#' @param meanX The average value of X. Defaults to \code{mean(X)}.
#'
#' @return A vector whose values are the kernel regression estimator.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_adjusted_est(X, 1:4, 1:3, 0.1, "gaussian", c(), FALSE, TRUE)
compute_adjusted_est <- function(X, x, t, h, kernel_name="gaussian", kernel_params=c(), custom_kernel = FALSE, pd = TRUE, type='covariance', meanX = mean(X)) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(x), length(x) >= 1, !any(is.na(x)),
            length(meanX) == 1, is.numeric(meanX), !is.na(meanX), !any(is.na(t)), is.numeric(t),
            length(t) >= 1, length(h) == 1, is.numeric(h), h > 0, is.logical(custom_kernel),
            type %in% c('covariance', 'correlation'), is.logical(pd))

  outer_x_x <- outer(x, x, '-')
  xij_mat <- Xij_mat(X)
  cov_vals <- c()
  if(custom_kernel) {
    for(tt in t) {
      K_ij <- get(kernel_name)(tt - outer_x_x, h, kernel_params)
      cov_vals <- c(cov_vals, sum(K_ij * xij_mat) / sum(K_ij))
    }
  }
  else {
    stopifnot(kernel_name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    for(tt in t) {
      K_ij <- kernel_symm(tt - outer_x_x, kernel_name, c(h, kernel_params[1], kernel_params[2]))
      cov_vals <- c(cov_vals, sum(K_ij * xij_mat) / sum(K_ij))
    }
  }

  if(pd) {
    # Perform DCT and set frequencies whose associated value is <0 to 0.
    cov_dct <- dct_1d(cov_vals)
    cov_dct[which(cov_dct < 0)] <- 0

    # Invert and return.
    cov_idct <- idct_1d(cov_dct)

    if(type == 'correlation') {
      cov_idct <- cov_idct / cov_idct[1]
    }

    return(cov_idct)
  }

  if(type == 'correlation') {
    cov_vals <- cov_vals / cov_vals[1]
  }
  return(cov_vals)
}

#' Make any function positive-definite
#'
#' This helper function can make any function positive-definite using methods proposed by P. Hall and his coauthors.
#'
#' @details
#' This function perform positive-definite adjustments proposed by P. Hall and his coauthors.
#'
#' Method 1 is as follows:
#' 1. Take the discrete cosine transform,
#' \eqn{\mathcal{F}^{c}(\hat{\rho}(t))}.
#' 2. Set all negative values to zero,
#' \eqn{\hat{\mathcal{F}}^{c}(\hat{\rho}(t)) = \mathcal{F}(\hat{\rho}(t))} if \eqn{\mathcal{F}(\hat{\rho}(t)) > 0} and \eqn{\hat{\mathcal{F}}^{c}(\hat{\rho}(t)) = 0} otherwise,
#' for all sample frequencies.
#' 3. Perform the inversion to obtain a new estimator.
#'
#' Method 2 is as follows:
#' 1. Take the discrete cosine transform
#' \eqn{\mathcal{F}^{c}(\hat{\rho}_{1}(t))}.
#' 2. Find the smallest frequency where its associated value in the spectral domain is negative
#' \deqn{\hat{\theta} = \inf\{ \theta > 0 :  \mathcal{F}^{c}(\hat{\rho}_{1}(t)) < 0\}.}
#' 3. Set all values starting at the frequency to zero.
#' 4. Perform the inversion.
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#'
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the nonparametric estimation of covariance functions. The Annals of Statistics (Vol. 22, Issue 4, pp. 2115–2134). 10.1214/aos/1176325774
#'
#' @param x A vector of numeric values - for example, an estimated autocovariance function.
#' @param method.1 Should method 1 be used (TRUE) or method 2 (FALSE).
#'
#' @return A vector that is the adjusted function.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' make_pd(X)
make_pd <- function(x, method.1 = TRUE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), is.logical(method.1))
  dct_x <- dct_1d(x)
  if(method.1) {
    dct_x[dct_x < 0] <- 0
  }
  else if(!method.1) {
    firstNeg <- which(dct_x[-1] < 0)[1] + 1

    if(length(firstNeg) == 0 || is.na(firstNeg)) {
      return(x)
    }

    # Set values after the first negative index onwards
    dct_x[firstNeg:length(dct_x)] <- 0
  }
  idct_x <- idct_1d(dct_x)
  return(idct_x)
}
