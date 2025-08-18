#' Generate Spline Knots.
#'
#' A helper function that generates \eqn{m + 2} spline knots of the form:
#' \deqn{
#' \kappa_{0} = 0 , \kappa_{1} = 1 / (m + 1) , \dots , \kappa_{m} = m / (m + 1) , \kappa_{m + 1} = 1 .
#' }
#' The knots are equally spaced with boundary knots \eqn{\kappa_{0} = 0} and \eqn{\kappa_{m + 1} = 1 .}
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param m The number of nonboundary knots.
#'
#' @return A numeric vector representing the knots, including the boundary knots.
#' @export
#'
#' @examples
#' generate_knots(3)
generate_knots <- function(m) {
  stopifnot(is.numeric(m), m > 0, m %% 1 == 0)
  knotVec <- c(0)
  for(i in 1:m) {
    knotVec <- c(knotVec, i / (m + 1))
  }
  knotVec <- c(knotVec, 1)
  return(knotVec)
}

#' Get a Specific \eqn{\tau_{i}}.
#'
#' A helper function that transforms the knots from [generate_knots] into the following form:
#' For \eqn{i = -p , -p + 1, \dots , -2, -1 , m + 2, m + 3, \dots , m + p , m + p + 1,} it is equal to \eqn{\tau_{i} = i / (m + 1)}, and for \eqn{i = 0, \dots , m + 1,} it is \eqn{\tau_{i} = \kappa_{i}.}
#' See Choi, Li & Wang (2013) page 615 for details.
#' This is a helper function of [get_taus].
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param i The knot index (\eqn{-p} through \eqn{m + p + 1}).
#' @param p The order of the splines.
#' @param m The number of nonboundary knots.
#' @param kVec Knot vector - see [generate_knots].
#'
#' @return The numerical value of \eqn{\tau_{i}.}
#' @export
#'
#' @examples
#' kVec <- generate_knots(2)
#' get_tau(1, 3, 2, kVec)
get_tau <- function(i, p, m, kVec) {
  stopifnot(is.numeric(p), p >= 0, p %% 1 == 0, is.numeric(m), m > 0, m %% 1 == 0,
            is.vector(kVec), length(kVec) > 0, all(!is.na(kVec)), all(kVec <= 1 & kVec >= 0))
  if(((i >= -p) & (i <= -1)) | (i >= (m + 2)) & (i <= (m + p + 1))) {
    return(i / (m + 1))
  }
  else if((i >= 0) & (i <= (m + 1))) {
    return(kVec[i + 1])
  }
  else {
    return(NA)
  }
}
#' Get all \eqn{\tau}.
#'
#' A helper function that repeatedly calls [get_tau] to obtain all \eqn{\tau_{i}, i=-p, \dots, m + p + 1,} where each \eqn{\tau_{i}} is as follows.
#' For \eqn{i = -p , -p + 1, \dots , -2, -1 , m + 2, m + 3, \dots , m + p , m + p + 1,} it is equal to \eqn{\tau_{i} = i / (m + 1)}, and for \eqn{i = 0, \dots , m + 1,} it is \eqn{\tau_{i} = \kappa_{i}.}
#' See Choi, Li & Wang (2013, p. 615) for details.
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param p The order of the splines.
#' @param m The number of nonboundary knots.
#'
#' @return A numeric vector of all \eqn{\tau_{i}, i = -p, \dots, m + p + 1.}
#' @export
#'
#' @examples
#' get_taus(3, 2)
get_taus <- function(p, m) {
  stopifnot(is.numeric(p), p >= 0, p %% 1 == 0, is.numeric(m), m > 0, m %% 1 == 0)
  kVec <- generate_knots(m)
  tauVec <- list()
  for(i in -p:(m + p + 1)) {
    tauVec[[paste0("", i, "")]] <- get_tau(i, p, m, kVec)
  }
  return(tauVec)
}

#' Compute Adjusted Splines.
#'
#' A helper function that is an implementation of the formula from Choi, Li & Wang (2013, p. 616),
#' \deqn{
#' f_{j}^{(l)}(x) = \frac{m + 1}{l} \left( f_{j}^{(l - 1)}(x + 1) - \tau_{j - p} f_{j}^{(l - 1)}(x) + \tau_{j - p + l + 1} f_{j + 1}^{(l - 1)}(x) -  f_{j + 1}^{(l - 1)}(x + 1)  \right) ,
#' }
#' where \eqn{m} is the number of nonboundary knots, \eqn{p} is the order of the spline, \eqn{l} is the order of the adjusted spline (the function \eqn{f_{j}^{(l)}(\cdot)}) and \eqn{j = 1, 2, \dots , m + p.}
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param x Argument of the function.
#' @param j Index of basis function of order \eqn{l}.
#' @param l Order of function.
#' @param p The order of the splines.
#' @param m The number of nonboundary knots.
#' @param taus Vector of \eqn{\tau}s, see [get_tau].
#'
#' @return A numeric value of the adjusted spline \eqn{f_{j}^{(l)}(x).}
#' @export
#'
#' @examples
#' taus <- get_taus(3, 2)
#' adjusted_spline(1, 2, 1, 3, 2, taus)
adjusted_spline <- function(x, j, l, p, m, taus) {
  stopifnot(is.numeric(x), j %% 1 == 0, is.numeric(j), j > 0, l %% 1 == 0, is.numeric(l), is.numeric(p), p >= 0, p %% 1 == 0,
            is.numeric(m), m > 0, m %% 1 == 0, is.vector(taus), length(taus) > 0, all(!is.na(taus)))
  # base case
  if(l == 0) {
    tau1 <- taus[[paste0("", j - p, "")]]
    tau2 <- taus[[paste0("", j - p + 1, "")]]
    if((tau1 >= 0 & tau1 <= 1) & (tau2 >= 0 & tau2 <= 1) ) {
      constant <- (m + 1) / (x + 1)
      inner <- tau2^(x + 1) - tau1^(x + 1)
      return(constant * inner)
    }
    return(0)
  }
  else {
    constant <- (m + 1) / l
    firstInner <- adjusted_spline(x + 1, j, l - 1, p, m, taus)
    tau_jp <- taus[[paste0("", j - p, "")]]
    secondInner <- tau_jp * adjusted_spline(x, j, l - 1, p, m, taus)

    tau_jpl1 <- taus[[paste0("", j - p + l + 1, "")]]

    thirdInner <- tau_jpl1 * adjusted_spline(x, j + 1, l - 1, p, m, taus)
    fourthInner <- adjusted_spline(x + 1, j + 1, l - 1, p, m, taus)

    retVal <- constant * (firstInner - secondInner + thirdInner - fourthInner)
    return(retVal)
  }
}

#' Construct Data Frame of Basis Functions.
#'
#' This helper function constructs a data frame with the following structure:
#' * One column for the x-values
#' * m + p columns values of squared basis functions evaluated at the corresponding x.
#'
#' @param x A vector of lags.
#' @param p The order of the splines.
#' @param m The number of nonboundary knots.
#' @param taus Vector of \eqn{\tau}s, see [get_tau].
#'
#' @return A data frame of the above structure.
#' @export
#'
#' @examples
#' taus <- get_taus(3, 2)
#' splines_df(seq(0, 2, by=0.25), 3, 2, taus)
splines_df <- function(x, p, m, taus) {
  stopifnot(is.numeric(x), length(x) > 0, all(!is.na(x)), is.numeric(p), p >= 0, p %% 1 == 0, is.numeric(m), m > 0, m %% 1 == 0,
            is.vector(taus), length(taus) > 0, all(!is.na(taus)))
  valsDf <- data.frame('lags' = x)
  for(j in 1:(m+p)) {
    jStr <- paste0('j', j)
    valsDf[, jStr] <- rep(NA, nrow(valsDf))
  }

  l <- p - 1
  for(i in 1:nrow(valsDf)) {
    for(j in 1:(m+p)) {
      valsDf[i, j + 1] <- adjusted_spline(valsDf[i, 1]^2, j, l, p, m, taus)
    }
  }

  return(valsDf)
}

#' Objective Function for WLS.
#'
#' This is the objective function to find the weights for each basis function in the minimising spline, see Choi, Li & Wang (2013, p. 617).
#' The parameters must be nonnegative, so a penalty of \eqn{10^{12}} is given if any parameters are negative.
#' The weights are chosen as per Choi, Li & Wang (2013, p. 617).
#'
#' @details
#' Let \eqn{\mathbf{\beta} = (\beta_{0}, \dots, \beta_{m + p})^{\prime}} be a vector of model coefficients, \eqn{\{f_{1}^{(p - 1)} , \dots , f_{m + p}^{(p - 1)} \}}
#' be a set of completely monotone basis functions, and \eqn{\widehat{C}(\cdot)} be an estimated covariance function.
#' As per Choi, Li & Wang (2013, p. 617), \eqn{\mathbf{\beta}} can be estimated via weighted-least squares,
#' \deqn{\
#' \hat{\mathbf{\beta}}_{WLS} = {\arg\min}_{\beta_{j} \ge 0} \sum_{i=1}^{L} w_{i} \left(\widehat{C}(h_{i}) - \sum_{j = 1}^{m + p} \beta_{j} f_{j}^{(p - 1)}(h_{i}^{2})  \right)^{2} ,
#' }
#' where \eqn{\{h_{1} , \dots , h_{L} \}} is a set of lags and \eqn{\{w_{1}, \dots , w_{L} \}} is a set of weights.
#' The set of weights is calculated in [splines_est], and they are of the form \eqn{w_{i} = (N - h_{i}) / ((1 - \widehat{C}(h_{i}))^{2}).}
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param par A vector of initial parameters to used in the minimisation process.
#' @param splines_df A data frame whose structure is defined in [splines_df].
#' @param weights A vector of weights, see the description.
#'
#' @return The value of the objective function at those parameters.
#' @export
#'
#' @examples
#' taus <- get_taus(3, 2)
#' x <- seq(0, 2, by=0.25)
#' maxLag <- 4
#' splines_df <- splines_df(x[1:maxLag], 3, 2, taus)
#' splines_df$estCov <- exp(-splines_df$lags^2) + 0.001
#' # pars are the inital parameters used in the minimisation process.
#' pars <- c(0.5, 0.5, 0.5, 0.5, 0.5)
#' weights <- c()
#' X <- rnorm(50)
#' for(i in 0:(maxLag - 1)) {
#' weights <- c(weights, (length(X) - i) / ( (1 - splines_df$estCov[i + 1])^2 ))
#' }
#' solve_spline(pars, splines_df, weights)
solve_spline <- function(par, splines_df, weights) {
  if(any(par < 0)) {
    return(10^12)
  }

  splineParVal <- 0
  for(i in 2:(ncol(splines_df) - 1)) {
    splineParVal <- splineParVal + par[i - 1] * splines_df[, i]
  }

  return( sum(weights * (splines_df$estCov - splineParVal)^2) )
}

#' Compute the Splines Estimator.
#'
#' Compute the estimated covariance function by using the method from Choi, Li & Wang (2013, pp. 614-617).
#' \deqn{C(\tau) = \sum_{j = 1}^{m + p} \beta_{j} f_{j}^{(p-1)}(\tau^{2}),}
#' where \eqn{m} is the number of nonboundary knots, \eqn{p} is the order of the splines, \eqn{\tau} is the isotropic distance, \eqn{\beta_{j}} are nonnegative weights and \eqn{f_{j}^{(p)}} are basis functions of order \eqn{p.}
#' For optimisation, the Nelder-Mead and L-BFGS-B methods are used, the one which selects parameters which minimises the objective function is chosen.
#'
#' @references Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611-630. 10.1007/s13253-013-0152-z
#'
#' @param X A vector representing observed values of the time series.
#' @param x A vector of lags.
#' @param estCov An estimated autocovariance function to fit to (a vector).
#' @param p The order of the splines.
#' @param m The number of nonboundary knots.
#' @param maxLag An optional parameter that determines the maximum lag to compute the estimated autocovariance function at. Defaults to \code{length(X) - 1}.
#' @param inital_pars An optional vector of parameters - can be used to fine tune the fit. By default, it is a vector of 0.5 whose length is \eqn{m+p.}
#' @param control An optional list of optimisation parameters used in the optimisation process, see \code{control} in [stats::optim].
#' @param type Compute either the 'autocovariance' or 'autocorrelation'. Defaults to 'autocovariance'.
#'
#' @return A vector whose values are the spline autocovariance estimates.
#' @export
#'
#' @importFrom stats optim
#' @importFrom stats runif
#'
#' @examples
#' X <- rnorm(100)
#' x <- seq(0, 5, by = 0.25)
#' maxLag <- 5
#' estCov <- standard_est(X, maxLag = maxLag)
#' estimated <- splines_est(X, x, estCov, 3, 2, maxLag = maxLag)
#' estimated
splines_est <- function(X, x, estCov, p, m, maxLag = length(X) - 1, type = "autocovariance", inital_pars = c(), control = list('maxit' = 1000)) {
  stopifnot(is.numeric(X), is.vector(X), all(!is.na(X)), is.numeric(x), is.vector(x), all(!is.na(x)),
            is.numeric(maxLag), maxLag >= 0, maxLag <= (length(estCov) - 1), is.numeric(estCov),
            is.vector(estCov), all(!is.na(estCov)), is.numeric(p), p >= 0, p %% 1 == 0,
            is.numeric(m), m > 0, m %% 1 == 0, type %in% c('autocovariance', 'autocorrelation'))

  if(is.vector(inital_pars)) {
    stopifnot(length(inital_pars) == p + m)
  }
  taus <- get_taus(p, m)

  splines_df <- splines_df(x[1:(maxLag + 1)], p, m, taus)
  splines_df[, 'estCov'] <- estCov

  weights <- c()
  for(i in 0:maxLag) {
    weights <- c(weights, (length(X) - i) / ( (1 - estCov[i + 1])^2 ))
  }

  pars <- rep(0.5, ncol(splines_df) - 2)
  if(length(inital_pars) != 0 && length(inital_pars) == ncol(splines_df - 2)) {
    pars <- inital_pars
  }
  # Add two optim results, one for Nelder-Mead and one for BFGS, mimicking what is done in the optimx package.
  optim_splines <- NA
  optim_splines_bfgs <- NA

  optim_splines <- optim(par=pars, fn = solve_spline, splines_df = splines_df, weights = weights, method = 'Nelder-Mead', control = control)
  converge <- optim_splines$convergence == 0

  optim_splines_bfgs <- optim(par=pars, fn = solve_spline, splines_df = splines_df, weights = weights, method = 'L-BFGS-B',
                         lower=rep(0, ncol(splines_df) - 2), control = control)
  converge_bfgs <- optim_splines_bfgs$convergence == 0

  optim_splines_betas <- c()
  # Check if both converge
  if(converge && converge_bfgs) {
    isNelder <- ifelse(optim_splines$value < optim_splines_bfgs$value, TRUE, FALSE)
    if(isNelder) {
      optim_splines_betas <- optim_splines$par
    }
    else {
      optim_splines_betas <- optim_splines_bfgs$par
    }
  }

  # if Nelder-Mead converges
  else if(converge && !converge_bfgs) {
    optim_splines_betas <- optim_splines$par
  }

  # if BFGS converges
  else if(!converge && converge_bfgs) {
    optim_splines_betas <- optim_splines_bfgs$par
  }

  # Neither converge
  else {
    # stop("Failed to converge after five tries.")
    stop("Failed to converge. Try new parameters or increasing the number of iterations.")
  }

  optim_vals <- c()
  for(i in 1:length(x[1:(maxLag + 1)])) {
    tempSum <- 0
    for(j in 2:(ncol(splines_df) - 1)) {
      tempSum <- tempSum + (optim_splines_betas[j - 1] * splines_df[i, j])
    }
    optim_vals <- c(optim_vals, tempSum)
  }

  if(length(optim_vals) != length(estCov)) {
    stop("Length of output differs from input")
  }

  if(type == 'autocorrelation') {
    optim_vals <- optim_vals / optim_vals[1]
  }
  return(optim_vals)
}
