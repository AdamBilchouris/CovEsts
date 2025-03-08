% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/kernel_regression_estimator.R
\name{truncated_point}
\alias{truncated_point}
\title{Compute \eqn{\rho(t)} for the truncated kernel regression estimator.}
\usage{
truncated_point(
  x,
  meanX,
  t,
  T1,
  T2,
  h,
  xij_mat,
  rhoT1,
  kernel_name = "gaussian",
  kernel_params = c(),
  custom_kernel = F
)
}
\arguments{
\item{x}{A vector of indices.}

\item{meanX}{The average value of X.}

\item{t}{The value at which the covariance function is calculated at.}

\item{T1}{The first truncation point, \eqn{T_{1} > 0.}}

\item{T2}{The second truncation point, \eqn{T_{2} > T_{1} > 0.}}

\item{h}{Bandwidth parameter.}

\item{xij_mat}{The matrix of pairwise covariance values.}

\item{rhoT1}{The value of the covariance function at T1.}

\item{kernel_name}{The name of the kernel function to be used. Possible values are:
"gaussian", "wave", "rational_quadratic", and "bessel_j". Alternatively, a custom kernel function can be provided, see \link{compute_corrected_standard_est}'s example.}

\item{kernel_params}{A vector of parameters of the kernel function. See \link{kernel_symm} for parameters.}

\item{custom_kernel}{If a custom kernel is to be used or not.}
}
\value{
The estimated covariance value at t.
}
\description{
This function computes \eqn{\rho(t)} for the truncated kernel regression estimator (see \link{compute_truncated_est}).
}
\examples{
X <- c(1, 2, 3, 4)
rhoT1 <- rho_T1(1:4, mean(X), 1, 0.1, Xij_mat(X),
                        "gaussian", c(), FALSE)
truncated_point(1:4, mean(X), 1, 0.1, 1, 0.1,
                        Xij_mat(X), rhoT1, "gaussian", c(), FALSE)
}
\references{
Hall, P., & Patil, P. (1994). Properties of nonparametric estimators of autocovariance for stationary random fields. In Probability Theory and Related Fields (Vol. 99, Issue 3, pp. 399–424). 10.1007/bf01199899

Hall, P., Fisher, N. I., & Hoffmann, B. (1994). On the Nonparametric Estimation of Covariance Functions. In The Annals of Statistics (Vol. 22, Issue 4). Institute of Mathematical Statistics. 10.1214/aos/1176325774
}
