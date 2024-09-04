compute_Xij_mat <- function(X) {
  stopifnot(length(X) >= 1, !any(is.na(X)))

  meanX <- mean(X)
  Xij <- matrix(nrow=length(X), ncol=length(X))
  for(i in 1:length(X)) {
    for(j in 1:length(X)) {
      Xij[i, ] <- (X[i] - meanX) * (X[j] - meanX)
    }
  }

  return(Xij)
}

compute_rho_T1 <- function(x, meanX, T1, h, Xij_mat, kernel="gaussian", kernel_params=c(), custom_kernel = F) {
  stopifnot(length(x) >= 1, !any(is.na(x)), is.double(meanX), T1 > 0, h > 0, is.matrix(Xij_mat), is.logical(custom_kernel))

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
    stopifnot(kernel %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    for(i in 1:length(x)) {
      tij <- x[i] - x
      t_tij <- T1 - tij
      X_ij <- Xij_mat[i, ]
      K_ij <- get(paste0("kernel_symm_", kernel))(t_tij, h, kernel_params[1], kernel_params[2])
      numerator <- K_ij * X_ij
      denominator <- K_ij

      numerators <- c(numerators, sum(numerator))
      denominators <- c(denominators, sum(denominator))
    }
  }

  return( ( sum(numerators) / sum(denominators) ) )
}

compute_truncated_point <- function(x, meanX, t, T1, T2, h, Xij_mat, rho_T1, kernel="gaussian",  kernel_params=c(), custom_kernel = F) {
  stopifnot(length(x) >= 1, !any(is.na(x)), is.double(meanX), !is.na(meanX), is.double(t), !is.na(t),
            is.double(T1), !is.nan(T1), T1 > 0, is.double(T2), !is.nan(T2), T2 > T1,
            h > 0, is.matrix(Xij_mat), !any(is.na(Xij_mat)), is.double(rho_T1), !is.na(rho_T1), is.logical(custom_kernel))
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
    stopifnot(kernel %in% c("gaussian", "wave", "rational_quadratic", "bessel_j"))
    # Case 1: 0 <= t <= T1
    # \hat{\rho}(t)
    if(t >= 0 && t <= T1) {
      for(i in 1:length(x)) {
        tij <- x[i] - x
        t_tij <- t - tij
        X_ij <- Xij_mat[i, ]
        K_ij <- get(paste0("kernel_symm_", kernel))(t_tij, h, kernel_params[1], kernel_params[2])
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

compute_truncated <- function(X, x, meanX, t, T1, T2, h, kernel="gaussian",  kernel_params=c(), custom_kernel = F) {
  stopifnot(length(x) >= 1, !any(is.na(x)), is.double(meanX), !is.na(meanX), is.double(t), !is.na(t),
            is.double(T1), !is.nan(T1), T1 > 0, is.double(T2), !is.nan(T2), T2 > T1,
            h > 0, is.logical(custom_kernel))

  Xij_mat <- compute_Xij_mat(X)
  rho_T1 <- compute_rho_T1(x, meanX, T1, h, Xij_mat, kernel, kernel_params, custom_kernel)

  vals_truncated_1 <- sapply(1:length(t), function(i) compute_truncated_point(x, meanX, t[i], T1, T2, h, Xij_mat, rho_T1, kernel, kernel_params, custom_kernel))
  return(vals_truncated_1)
}

