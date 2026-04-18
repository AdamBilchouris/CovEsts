#' Random Block Locations
#'
#' This function performs random sampling to obtain random starting locations for block bootstrap.
#'
#' @details
#' This function performs random sampling to obtain random starting locations for block bootstrap.
#' If \code{type = 'moving'}, the set \eqn{\{1, \dots, N - \ell  + 1\}} is randomly sampled, with replacement, \eqn{k} times to obtain random block locations for moving block bootstrap.
#' If \code{type = 'circular'}, the set \eqn{\{1, \dots, N\}} is randomly sampled, with replacement, \eqn{k} times to obtain random block locations for moving block bootstrap.
#'
#' @param N The length of the observation window.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @return A vector of length \code{k} whose values are random block locations.
#' @export
#'
#' @examples
#' starting_locs(4, 2, 2)
starting_locs <- function(N, l, k, boot_type = c("moving", "circular")) {
  stopifnot(is.numeric(N), N >= 1, !is.na(N), is.numeric(l), l >= 1, !is.na(l),
            is.numeric(k), k >= 1, !is.na(k))

  boot_type <- match.arg(boot_type)

  sampleVec <- switch(boot_type,
                      moving = 1:(N - l + 1),
                      circular = 1:N)

  return(sample(sampleVec, k, replace = TRUE))
}

#' Block Bootstrap Sample
#'
#' This function generates block bootstrap samples for either moving block bootstrap or circular bootstrap.
#'
#' @details
#' This function generates a block bootstrap sample for a time series \eqn{X.}
#' For the moving block bootstrap and circular bootstrap procedures see [block_bootstrap] and the included references.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @param X A vector representing observed values of the time series.
#' @param l The block length considered for bootstrap.
#' @param k The number of blocks considered for bootstrap.
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#'
#' @return A vector of length \code{length(X)} whose values are a bootstrapped time series.
#' @export
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' bootstrap_sample(X, 2, 3)
bootstrap_sample <- function(X, l, k, boot_type = c("moving", "circular")) {
  stopifnot(is.numeric(X), length(X) >= 1, !is.na(X), is.numeric(l), l >= 1, !is.na(l),
            is.numeric(k), k >= 1, !is.na(k))
  boot_type <- match.arg(boot_type)

  N <- length(X)
  S <- starting_locs(N, l, k, boot_type)
  block_offsets <- 0:(l-1)

  if(boot_type == 'moving') {
    indices <- outer(S, block_offsets, "+")
  }

  else if(boot_type == 'circular') {
    indices <- outer(S, block_offsets, "+")
    indices <- ((indices - 1) %% N) + 1
  }

  Y <- X[t(indices)]
  Y <- Y[1:N]
  return(Y)
}

#' Block Bootstrap
#'
#' This function performs block bootstrap (moving or circular) to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#'
#' @details
#' This function performs block bootstrap to obtain a \eqn{(1-\alpha)\%} confidence-interval for the autocovariance function. It will also compute average autocovariance function across all bootstrapped estimates.
#'
#' Moving block bootstrap can be described as follows.
#' For some times series \eqn{X(1), X(2), \dots, X(n),} construct \eqn{k \in N} overlapping blocks of the form
#' \eqn{B_{i} = (X(i), \dots, X(i + \ell - 1)),} where \eqn{\ell \in \{1, \dots , n\}} is the block length.
#' Randomly sample, with replacement, from the discrete uniform distribution with on \eqn{\{1, \dots, n - \ell + 1\}} to obtain a set of random starting locations \eqn{I_{1}, \dots, I_{k}.}
#' Construct a bootstrapped time series \eqn{B_{1}^{\ast}, B_{2}^{\ast}, \dots, B_{k}^{\ast},} where \eqn{B_{i}^{\ast} = B_{I_{i}}.}
#' The bootstrapped time series is truncated to have length \eqn{n,} and will be of the form \eqn{X^{\ast}(1), \dots , X^{\ast}(n).}
#' The autocovariance function is then computed for the bootstrapped time series.
#'
#' An alternative to moving block bootstrap is circular block bootstrap.
#' Circular block bootstrap uses the time series like a circle, that is, the observation at \eqn{n + i} is the same as the observation at location \eqn{i.}
#' For example, for the block \eqn{B_{n - \ell + 2}}, we obtain \eqn{(X(n - \ell + 2) , \dots , X(n), X(n + 1))} is the same as \eqn{(X(n - \ell + 2) , \dots , X(n), X(1)).}
#' When performing random sampling to obtain starting locations, the set \eqn{\{1, \dots, n\}} is now considered.
#' The procedure for constructing the bootstrap time series after constructing the blocks and selecting the starting indices is the same as moving block bootstrap.
#'
#' This process is repeated \code{n_bootstrap} times to obtain \code{n_boostrap} estimates of the autocovariance function using the bootstrapped time series, for which the average autocovariance function
#' and the \eqn{(1 - \alpha)\%} confidence intervals are constructed pointwise for each lag.
#'
#' The choice of the block length, \eqn{\ell ,} depends on the degree of dependence present in the time series. If the time series has a high degree of dependence, a larger block size should be chosen to ensure the dependency structure is maintained within the block.
#'
#' Any estimator of the autocovariance function can be used in this function, including a custom estimator not in the package, see the examples.
#'
#' When \code{parallel = TRUE}, L'Ecuyer-CMRG is used to ensure independent random numbers across workers. The serial version uses the Mersenne-Twister algorithm.
#'
#' @references
#' Chapters 2.5 and 2.7 in Lahiri, S. N. (2003). Resampling Methods for Dependent Data. Springer. https://doi.org/10.1007/978-1-4757-3803-2
#'
#' Künsch, H. R. (1989). The Jackknife and the Bootstrap for General Stationary Observations. The Annals of Statistics 17(3), 1217-1241. https://doi.org/10.1214/aos/1176347265
#'
#' Politis, D. N. & Romano, J. P. (1991). A Circular Block-Resampling Procedure for Stationary Data. In R. LePage & L. Billard, eds, Exploring the Limits of Bootstrap, Wiley, 263-270.
#'
#' @param X A vector representing observed values of the time series.
#' @param maxLag The maximum lag to compute the estimated autocovariance function at.
#' @param x A vector of lag indices. Defaults to the sequence \code{0:length(X)}. Must be at least as large as \code{maxLag + 1}.
#' @param n_bootstrap The number of times to run moving block bootstrap. Defaults to 100.
#' @param l The block length considered for bootstrap. Defaults to \eqn{\lceil N \rceil^{1/3}}, where \eqn{N} is the length of the observation window.
#' @param estimator The function name of the estimator to use. Defaults to \code{standard\_est}.
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#' @param alpha The quantile used to compute the \eqn{(1 - \alpha)\%} confidence interval. Defaults to \eqn{0.05.}
#' @param boot_type What type of block bootstrap should be used, either 'moving' for moving block bootstrap or 'circular' for circular block bootstrap.
#' @param parallel Whether or not the bootstrap computations should be done in parallel or not. Defaults to \code{FALSE}.
#' @param ncores The number of cores to be used in the parallel bootstrap computations. Defaults to the number cores - 1 (threads if hyperthreading is available), calculated from \code{parallel::detectCores() - 1}.
#' @param cl_export A vector of any additional functions or variables to export for parallel computations. This may be required if \code{estimator} is not within the package. Defaults to \code{NULL}.
#' @param boot_seed An integer seed for reproducibility. This is used for
#' @param cl An optional cluster object created by \code{parallel::makeCluster}. Defaults to \code{NULL}, which creates a temporary PSOCK cluster.
#' @param ... Optional named arguments to the chosen estimator. See the examples.
#'
#' @return A \code{BootEsts} S3 object (list) with the following values
#' \describe{
#'  \item{\code{acf_avg}}{A numeric vector containing the average autocovariance/autocorrelation bootstrap estimate.}
#'  \item{\code{lags}}{A numeric vector containing the lag indices used to compute the estimates on.}
#'  \item{\code{acf_orig}}{A numeric vector containing the nonbootstrapped autocovariance/autocorrelation estimate.}
#'  \item{\code{acf_mat}}{A numeric matrix of the a matrix of all of the bootstrap estimates.}
#'  \item{\code{conf_lower}}{A numeric vector containing the lower bounds for the estimated pointwise confidence interval.}
#'  \item{\code{conf_upper}}{A numeric vector containing the upper bounds for the estimated pointwise confidence interval.}
#'  \item{\code{est_type}}{The type of estimate, namely 'autocorrelation' or 'autocovariance', this depends on the argument \code{type}.}
#'  \item{\code{est_used}}{The estimator function used, this depends on the argument \code{estimator}.}
#'  \item{\code{boot_type}}{The type of estimate, namely 'moving' or 'circular', this depends on the argument \code{boot_type}.}
#'  \item{\code{alpha}}{A numeric value, which is the \eqn{\alpha} value used to compute the confidence intervals, this depends on the argument \code{alpha}.}
#' }
#' @export
#'
#' @importFrom stats quantile
#'
#' @examples
#' X <- c(1, 2, 3, 3, 2, 1)
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, type = 'autocorrelation')
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, type = 'autocovariance')
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, estimator=tapered_est,
#'     rho = 0.5, window_name = 'blackman', window_params = c(0.16),
#'     type='autocorrelation')
#
#' my_cov_est <- function(X, maxLag) {
#'   n <- length(X)
#'   covVals <- rep(0, maxLag + 1)
#'   for(h in 0:maxLag) {
#'     covVals_t <- (X[1:(n-h)] - mean(X)) * (X[(1+h):n] - mean(X))
#'     covVals[h] <- sum(covVals_t) / (n - h)
#'   }
#'   return(covVals)
#' }
#' block_bootstrap(X, 4, n_bootstrap = 3, l = 2, estimator=my_cov_est)
#'
#' plot(LakeHuron, main="Lake Huron Levels", ylab="Feet")
#' X <- as.vector(LakeHuron)
#' block_bootstrap(X, 20, n_bootstrap = 100, l = 40, type = 'autocorrelation')
#' block_bootstrap(X, 20, n_bootstrap = 100, l = 40, type = 'autocorrelation')
#' block_bootstrap(X, 20, n_bootstrap = 100, l = 40, estimator=tapered_est,
#'     rho = 0.5, window_name = 'blackman', window_params = c(0.16),
#'     type='autocorrelation')
#'
#' my_cov_est <- function(X, maxLag) {
#'   n <- length(X)
#'   covVals <- rep(0, maxLag + 1)
#'   for(h in 0:maxLag) {
#'     covVals_t <- (X[1:(n-h)] - mean(X)) * (X[(1+h):n] - mean(X))
#'     covVals[h] <- sum(covVals_t) / (n - h)
#'   }
#'   return(covVals)
#' }
#' block_bootstrap(X, 20, n_bootstrap = 100, l = 40, estimator = my_cov_est,
#'     type = 'autocorrelation')
#' \dontrun{
#' library(parallel)
#' X <- c(1, 2, 3, 3, 2, 1)
#' my_cl <- makePSOCKcluster(2)
#' block_bootstrap(X, 4, n_bootstrap = 1000, l = 3, parallel = TRUE, cl = my_cl)
#' stopCluster(my_cl)
#' }
block_bootstrap <- function(X, maxLag, x = 0:length(X), n_bootstrap = 100, l = ceiling(length(X)^(1/3)), estimator = standard_est, type = c("autocovariance", "autocorrelation"), alpha = 0.05, boot_type = c("moving", "circular"), parallel = FALSE, ncores = parallel::detectCores() - 1, cl_export = NULL, boot_seed = NULL, cl = NULL, ...) {
  stopifnot(is.numeric(X), length(X) >= 1, !any(is.na(X)), is.numeric(maxLag), length(maxLag) == 1,
            maxLag > 0, maxLag <= (length(X) - 1), maxLag %% 1 == 0, is.numeric(x), length(x) >= maxLag,
            is.numeric(n_bootstrap), n_bootstrap > 0, n_bootstrap %% 1 ==0, is.numeric(l), length(l) == 1,
            l > 0, l <= length(X), l %% 1 == 0, is.function(estimator), is.numeric(alpha), alpha <= 1, alpha >= 0,
            is.logical(parallel))

  type <- match.arg(type)
  boot_type <- match.arg(boot_type)

  ellipsis <- list(...)

  N <- length(X)
  k <- ceiling(N / l)
  acf_mat <- matrix(NA, ncol = maxLag + 1, nrow = n_bootstrap)

  boot_results <- NULL
  boot_seed <- if(is.null(boot_seed)) sample.int(.Machine$integer.max, 1) else boot_seed

  boot_iter <- function(i) {
    Y <- bootstrap_sample(X, l, k, boot_type)
    attempts <- 1
    max_attempts <- 100

    while(length(unique(Y)) == 1 && attempts < max_attempts) {
      Y <- bootstrap_sample(X, l, k, boot_type)
      attempts <- attempts + 1
    }

    if(attempts == max_attempts) {
      stop("Failed to generate a nonconstant bootstrap sample after 100 attempts. Exiting.")
    }

    est_args <- c(list(Y, maxLag = maxLag), ellipsis)
    return(as.numeric(do.call(estimator, est_args)))
  }

  if(parallel) {
    # See if the user objects actually exist before wasting time creating the PSOCK cluster.
    if(!is.null(cl_export)) {
      missing_exports <- cl_export[!sapply(cl_export, exists, envir = parent.frame())]
      if(length(missing_exports) > 0) {
        stop("The following objects in `cl_export` were not found in your environment: ", paste0(missing_exports, collapse = ", "))
      }
    }

    if(is.null(cl)) {
      cl <- parallel::makePSOCKcluster(ncores)
      on.exit(parallel::stopCluster(cl), add = TRUE)
    }

    varList <- c("X", "l", "k", "boot_type", "maxLag", "estimator", "ellipsis", "boot_iter")
    parallel::clusterExport(cl, varlist = varList, envir = environment())

    package_funcs <- c("bootstrap_sample", "starting_locs")
    parallel::clusterExport(cl, varlist = package_funcs, envir = asNamespace("CovEsts"))

    parallel::clusterSetRNGStream(cl, boot_seed)

    # Export user supplied objects
    parallel::clusterExport(cl, varlist = cl_export, envir = parent.frame())

    boot_results <- parallel::parLapply(cl, 1:n_bootstrap, boot_iter)
  }

  else {
    prev_rng <- RNGkind()
    on.exit(RNGkind(kind = prev_rng[1], normal.kind = prev_rng[2]), add = TRUE)
    set.seed(boot_seed)

    boot_results <- lapply(1:n_bootstrap, boot_iter)
  }

  acf_mat <- do.call(rbind, boot_results)

  original_acf <- as.numeric(estimator(X, ..., maxLag = maxLag))
  if(type == 'autocorrelation') {
    original_acf <- original_acf / original_acf[1]

    acf_mat <- acf_mat / acf_mat[, 1]
  }

  quantiles <- apply(acf_mat, 2, quantile, probs = c(alpha / 2, 1 - (alpha / 2)))
  lower_vals <- quantiles[1, ]
  upper_vals <- quantiles[2, ]
  average_acf <- colMeans(acf_mat)

  res <- list(acf_avg = average_acf, lags=x[1:(maxLag + 1)], acf_orig = original_acf,
              acf_mat = acf_mat, conf_lower = lower_vals, conf_upper = upper_vals,
              est_type = type, est_used = paste0(substitute(estimator)),
              boot_type = boot_type, alpha = alpha)

  return(structure(res, class = "BootEsts"))
}
