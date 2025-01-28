#' Compute \eqn{X_{ij}} Matrix
#'
#' This computes the matrix of pairwise covariance values, \eqn{X_{ij},} for the kernel regression estimator,
#' \deqn{X_{ij} = (X_{i} - \bar{X}) (X_{j} - \bar{X}) .}
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. In Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the Nonparametric Estimation of Covariance Functions. In The Annals of Statistics (Vol. 22, Issue 4). Institute of Mathematical Statistics. 10.1214/aos/1176325774
#'
#' @param X A vector of values
#'
#' @return A matrix of size \eqn{N \times N}, where \eqn{N} is the length of the vector X.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_Xij_mat(X)
compute_Xij_mat <- function(X) {
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
#' This function computes the truncated kernel regression estimator at \eqn{\rho(T_{1}).}
#' This saves repeated computations in compute_truncated.
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. In Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the Nonparametric Estimation of Covariance Functions. In The Annals of Statistics (Vol. 22, Issue 4). Institute of Mathematical Statistics. 10.1214/aos/1176325774
#'
#' @param x A vector of indices.
#' @param meanX The mean of the process \eqn{X.}
#' @param T1 The first trunctation point.
#' @param h Bandwidth parameter.
#' @param Xij_mat The matrix of pairwise covariance values.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j".
#' @param kernel_params Any parameters for the kernel function.
#' @param custom_kernel A boolean determining whether or not a custom kernel is used.
#'
#' @return The empirical covariance function at \eqn{T_{1}}.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_rho_T1(1:4, mean(X), 1, 0.1, compute_Xij_mat(X), "gaussian", c(), FALSE)
compute_rho_T1 <- function(x, meanX, T1, h, Xij_mat, kernel_name="gaussian", kernel_params=c(), custom_kernel = FALSE) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            length(T1) == 1, is.numeric(T1), !is.na(T1), T1 > 0, length(h) == 1, is.numeric(h), !is.na(h), h > 0,
            is.numeric(Xij_mat), is.matrix(Xij_mat), !any(is.na(Xij_mat)), is.logical(custom_kernel))

  numerators <- c()
  denominators <- c()

  # Duplicating code as to only have 1 if statement, not length(x) if statements.
  if(custom_kernel) {
    for(i in 1:length(x)) {
      tij <- x[i] - x
      t_tij <- T1 - tij
      X_ij <- Xij_mat[i, ]
      K_ij <- get(kernel)(t_tij, h, kernel_params)
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
      X_ij <- Xij_mat[i, ]
      # K_ij <- get(paste0("kernel_symm_", kernel))(t_tij, h, kernel_params[1], kernel_params[2])
      K_ij <- get("kernel_symm")(t_tij, kernel_name, c(h, kernel_params[1], kernel_params[2]))
      numerator <- K_ij * X_ij
      denominator <- K_ij

      numerators <- c(numerators, sum(numerator))
      denominators <- c(denominators, sum(denominator))
    }
  }

  return( ( sum(numerators) / sum(denominators) ) )
}

#' Compute \eqn{\rho(t)} for the truncated kernel regression estimator.
#'
#' This function computes the truncated kernel regression estimator at \eqn{\rho(t).}
#'
#' @references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. In Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the Nonparametric Estimation of Covariance Functions. In The Annals of Statistics (Vol. 22, Issue 4). Institute of Mathematical Statistics. 10.1214/aos/1176325774
#'
#' @param x A vector of indices.
#' @param meanX The mean of the process \eqn{X.}
#' @param t The value at which the covariance function is calculated at.
#' @param T1 The first truncation point, \eqn{T_{1} > 0.}
#' @param T2 The second truncation point, \eqn{T_{2} > T_{1} > 0.}
#' @param h Bandwidth parameter.
#' @param Xij_mat The matrix of pairwise covariance values.
#' @param rho_T1 The value of the covariance function at T1.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j".
#' @param kernel_params Any parameters for the kernel function.
#' @param custom_kernel A boolean determining whether or not a custom kernel is used.
#'
#' @return The covariance value at t.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' rhoT1 <- compute_rho_T1(1:4, mean(X), 1, 0.1, compute_Xij_mat(X),
#'                         "gaussian", c(), FALSE)
#' compute_truncated_point(1:4, mean(X), 1, 0.1, 1, 0.1,
#'                         compute_Xij_mat(X), rhoT1, "gaussian", c(), FALSE)
compute_truncated_point <- function(x, meanX, t, T1, T2, h, Xij_mat, rho_T1, kernel_name="gaussian",  kernel_params=c(), custom_kernel = F) {
  stopifnot(is.numeric(x), length(x) >= 1, !any(is.na(x)), length(meanX) == 1, is.numeric(meanX), !is.na(meanX),
            length(t) == 1, is.numeric(t), !is.na(t), length(T1) == 1, is.numeric(T1), !is.na(T1), T1 > 0,
            length(T2) == 1, is.numeric(T2), !is.na(T2), T2 > T1, length(h) == 1, is.numeric(h), h > 0, !is.na(h),
            is.numeric(Xij_mat), is.matrix(Xij_mat), !any(is.na(Xij_mat)), length(rho_T1) == 1, is.numeric(rho_T1),
            !is.na(rho_T1), is.logical(custom_kernel))
  numerators <- c()
  denominators <- c()


  if(custom_kernel) {
    # Case 1: 0 <= t <= T1
    # \hat{\rho}(t)
    if(t >= 0 && t <= T1) {
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t - tij
        X_ij <- Xij_mat[i, ]
        K_ij <- get(kernel)(t_tij, h, kernel_params)
        numerator <- K_ij * X_ij
        denominator <- K_ij

        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      return( sum(numerators) / sum(denominators) )
    }

    # Case 2: T1 < t <= T2
    # \hat{\rho}(T1) (T2 - t) / (T2 - T1)
    else if(T1 < t && t <= T2) {
      linear_part <- (T2 - t) * (T2 - T1)^(-1)
      return( rho_T1 * linear_part )
    }

    else {
      return(0)
    }
  }

  else {
    stopifnot(kernel_name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    # Case 1: 0 <= t <= T1
    # \hat{\rho}(t)
    if(t >= 0 && t <= T1) {
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t - tij
        X_ij <- Xij_mat[i, ]
        # K_ij <- get(paste0("kernel_symm_", kernel))(t_tij, h, kernel_params[1], kernel_params[2])
        K_ij <- get("kernel_symm")(t_tij, kernel_name, c(h, kernel_params[1], kernel_params[2]))
        numerator <- K_ij * X_ij
        denominator <- K_ij

        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      return( sum(numerators) / sum(denominators) )
    }

    # Case 2: T1 < t <= T2
    # \hat{\rho}(T1) (T2 - t) / (T2 - T1)
    else if(T1 < t && t <= T2) {
      linear_part <- (T2 - t) * (T2 - T1)^(-1)
      return( rho_T1 * linear_part )
    }

    else {
      return(0)
    }
  }
}

#' Compute 1D Discrete Cosine Transform
#'
#' This computes the Type-II discrete cosine transform.
#'
#' @param X A vector of values for which the discrete cosine transform is being computed.
#'
#' @return A vector whose values are the transform.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3)
#' compute_1d_dct(X)
compute_1d_dct <- function(X) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)))
  newX <- c(X, rev(X))
  zerosX <- rep(0, length(newX))
  newX <- c(rbind(zerosX, newX))

  dct_seq <- 0.5*Re(stats::fft(newX))[1:(length(newX)/4)]
  return(dct_seq)
}

#' Compute 1D Inverse Discrete Cosine Transform
#'
#' This computes the inverse of the Type-II discrete cosine transform.
#'
#' @param X A vector of values for which the discrete cosine transform is being computed.
#'
#' @return A vector whose values are the inverse transform.
#' @export
#'
#' @examples
#' X <- compute_1d_dct(c(1, 2, 3))
#' compute_1d_idct(X)
compute_1d_idct <- function(X) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)))
  # Reconstruct original frequency (i.e. the full dct_seq from compute_1d_dct)
  dct_full <- c(X, 0, -rev(X), -X[-1], 0, rev(X[-1]))

  # Perform the inversion
  idct_vals <- Re(stats::fft(dct_full, inverse = T)) * (2 / length(dct_full))

  # Recover original signal, drop 0s
  idct_vals_2 <- idct_vals[1:(2*length(X))]
  idct_vals_2 <- idct_vals_2[idct_vals_2 != 0]

  return(idct_vals_2)
}

#' Compute the Truncated Kernel Regression Estimator.
#'
#' \deqn{\hat{\rho}_{1}(t) = \begin{cases}
#' \hat{\rho}(t) & 0 \leq t \leq T_{1} \\
#' \hat{\rho}(T_{1}) (T_{2} - t)(T_{2} - T_{1})^{-1} & T_{1} < t \leq T_{2} \\
#' 0 & t > T_{2}
#'\end{cases}}
#'where \eqn{\hat{\rho}(\cdot)} is the kernel regression estimator.
#'
#' When making the estimator positive-definite, the following is done:
#' 1. Take the discrete cosine transform.
#' 2. Find the smallest frequency where its associated value in the spectral domain is negative.
#' 3. Set all values starting at the frequency to zero.
#' 4. Perform the inversion.
#'
#' This ensures the covariance function estimate is positive-definite.
#'
#'@references
#' Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. In Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899
#' Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the Nonparametric Estimation of Covariance Functions. In The Annals of Statistics (Vol. 22, Issue 4). Institute of Mathematical Statistics. 10.1214/aos/1176325774
#'
#' @param X A vector of values representing a stochastic process.
#' @param x A vector of indices.
#' @param meanX The mean of the process \eqn{X.}
#' @param t The values at which the covariance function is calculated at.
#' @param T1 The first truncation point, \eqn{T_{1} > 0.}
#' @param T2 The second truncation point, \eqn{T_{2} > T_{1} > 0.}
#' @param h Bandwidth parameter.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j".
#' @param kernel_params Any parameters for the kernel function.
#' @param custom_kernel A boolean determining whether or not a custom kernel is used.
#' @param pd Whether or not the estimator returned is positive-definite.
#'
#' @return A vector whose values are the estimated truncated kernel regression estimator.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_truncated(X, 1:4, mean(X), 1:3, 1, 2, 0.1,
#'                   "gaussian", c(), FALSE, TRUE)
compute_truncated <- function(X, x, meanX, t, T1, T2, h, kernel_name="gaussian",  kernel_params=c(), custom_kernel = FALSE, pd = TRUE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), length(x) >= 1, is.numeric(x), !any(is.na(x)),
            length(meanX) == 1, is.numeric(meanX), !is.na(meanX), is.numeric(t), length(t) >= 1,
            length(T1) == 1, is.numeric(T1), !is.na(T1), T1 > 0, length(T2) == 1, is.numeric(T2),
            !is.na(T2), T2 > T1, length(h) == 1, is.numeric(h), h > 0,
            is.logical(custom_kernel), is.logical(pd))

  Xij_mat <- compute_Xij_mat(X)
  rho_T1 <- compute_rho_T1(x, meanX, T1, h, Xij_mat, kernel_name, kernel_params, custom_kernel)

  # vals_truncated_1 <- sapply(1:length(t), function(i) compute_truncated_point(x, meanX, t[i], T1, T2, h, Xij_mat, rho_T1, kernel, kernel_params, custom_kernel))
  vals_truncated_1 <- c()
  for(i in 1:length(t)) {
    vals_truncated_1 <- c(vals_truncated_1, compute_truncated_point(x, meanX, t[i], T1, T2, h, Xij_mat, rho_T1, kernel_name, kernel_params, custom_kernel))
  }

  if(pd) {
    # Perform DCT
    vals_truncated_1_dct <- compute_1d_dct(vals_truncated_1)

    # Find the first nonzero negative index.
    firstNeg <- which(vals_truncated_1_dct[-1] < 0)[1] + 1

    if(length(firstNeg) == 0 || is.na(firstNeg)) {
      return(vals_truncated_1)
    }

    # Set values after the first negative index onwards
    vals_truncated_1_dct[firstNeg:length(vals_truncated_1_dct)] <- 0

    # Inversion
    vals_truncated_1_idct <- compute_1d_idct(vals_truncated_1_dct)

    return(vals_truncated_1_idct)
  }

  return(vals_truncated_1)
}

#' Compute the Kernel Regression Estimator.
#'
#' \deqn{
#' \hat{\rho}(t) = \left( \sum_{i} \sum_{j} \check{X}_{ij} K((t - t_{ij}) / h) \right) \left( \sum_{i} \sum_{j} K((t - t_{ij}) / h) \right)^{-1},
#' }
#' where \eqn{\check{X}_{ij} = (X(t_{i}) - \bar{X}) (X(t_{j}) - \bar{X})}, and \eqn{t_{ij} = t_{i} - t_{j}.}
#'
#' When making the estimator positive-definite, the following is done:
#' 1. Take the discrete cosine transform.
#' 2. Set all negative values to zero.
#' 3. Perform the inversion.
#' This ensures the covariance estimate is positive-definite.
#'
#' @param X A vector of values representing a stochastic process.
#' @param x A vector of indices.
#' @param meanX The mean of the process \eqn{X.}
#' @param t The values at which the covariance function is calculated at.
#' @param h Bandwidth parameter.
#' @param kernel_name The name of the kernel function to be used. Possible values are:
#' "gaussian", "wave", "rational_quadratic", and "bessel_j".
#' @param kernel_params Any parameters for the kernel function.
#' @param custom_kernel A boolean determining whether or not a custom kernel is used.
#' @param pd Whether or not the estimator returned is positive-definite.
#'
#' @return A vector whose values are the estimated kernel regression estimator.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 4)
#' compute_adjusted(X, 1:4, mean(X), 1:3, 0.1, "gaussian", c(), FALSE, TRUE)
compute_adjusted <- function(X, x, meanX, t, h, kernel_name="gaussian", kernel_params=c(), custom_kernel = FALSE, pd = TRUE) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(x), length(x) >= 1, !any(is.na(x)),
            length(meanX) == 1, is.numeric(meanX), !is.na(meanX), !any(is.na(t)), is.numeric(t),
            length(t) >= 1, length(h) == 1, is.numeric(h), h > 0, is.logical(custom_kernel), is.logical(pd))

  Xij_mat <- compute_Xij_mat(X)
  cov_vals <- c()
  if(custom_kernel) {
    for(ti in 1:length(t)) {
      numerators <- c()
      denominators <- c()
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t[ti] - tij
        X_ij <- Xij_mat[i, ]
        K_ij <- get(kernel_name)(t_tij, h, kernel_params)
        numerator <- K_ij * X_ij
        denominator <- K_ij

        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      cov_vals <- c(cov_vals, sum(numerators) / sum(denominators))
    }
  }
  else {
    stopifnot(kernel_name %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    for(ti in 1:length(t)) {
      numerators <- c()
      denominators <- c()
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t[ti] - tij
        X_ij <- Xij_mat[i, ]
        # K_ij <- get(paste0("kernel_symm_", kernel))(t_tij, h, kernel_params[1], kernel_params[2])
        K_ij <- get("kernel_symm")(t_tij, kernel_name, c(h, kernel_params[1], kernel_params[2]))
        numerator <- K_ij * X_ij
        denominator <- K_ij

        numerators <- c(numerators, sum(numerator))
        denominators <- c(denominators, sum(denominator))
      }
      cov_vals <- c(cov_vals, sum(numerators) / sum(denominators))
    }
  }

  if(pd) {
    # Perform DCT and set frequencies whose associated value is <0 to 0.
    cov_dct <- compute_1d_dct(cov_vals)
    cov_dct[which(cov_dct < 0)] <- 0

    # Invert and return.
    cov_idct <- compute_1d_idct(cov_dct)

    return(cov_idct)
  }
  return(cov_vals)
}
