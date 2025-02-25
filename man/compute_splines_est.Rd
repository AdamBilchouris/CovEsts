% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/B_splines_est.R
\name{compute_splines_est}
\alias{compute_splines_est}
\title{Compute the Splines Estimator.}
\usage{
compute_splines_est(
  X,
  x,
  maxLag,
  estCov,
  p,
  m,
  inital_pars = c(),
  control = list(maxit = 1000)
)
}
\arguments{
\item{X}{A vector representing the process.}

\item{x}{A vector of lags.}

\item{maxLag}{The maximum lag to compute the autocovariance function at.}

\item{estCov}{An estimated covariance function to fit to (a vector).}

\item{p}{The order of the splines.}

\item{m}{The number of nonboundary knots.}

\item{inital_pars}{An optional vector of parameters - can be used to fine tune the fit. By default, it is a vector of 0.5 whose length is \eqn{m+p.}}

\item{control}{An optional list of optimisation parameters used in the optimisation process, see \code{control} in \link[stats:optim]{stats::optim}.}
}
\value{
A vector whose values are the spline autocovariance estimator.
}
\description{
Compute the estimated covariance function found in Choi, Li & Wang (2013, pp. 614-617).
\deqn{C(\tau) = \sum_{j = 1}^{m + p} \beta_{j} f_{j}^{[p]}(\tau^{2}),}
where \eqn{m} is the number of nonboundary knots, \eqn{p} is the order of the splines, \eqn{\tau} is the isotropic distance, \eqn{\beta_{j}} are nonnegative weights and \eqn{f_{j}^{[p]}} are basis functions of order \eqn{p.}
For optimisation, the Nelder-Mead and L-BFGS-B methods are used, the one which selects parameters which minimises the objective function is chosen.
}
\examples{
X <- rnorm(100)
x <- seq(0, 5, by = 0.25)
maxLag <- 5
estCov <- compute_standard_est(X, maxLag)
estimated <- compute_splines_est(X, x, maxLag, estCov, 3, 2)
estimated
}
\references{
Choi, I., Li, B. & Wang, X. Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES 18, 611–630 (2013). https://doi.org/10.1007/s13253-013-0152-z
}
