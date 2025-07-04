% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/B_splines_est.R
\name{generate_knots}
\alias{generate_knots}
\title{Generate Spline Knots.}
\usage{
generate_knots(m)
}
\arguments{
\item{m}{The number of nonboundary knots.}
}
\value{
A numeric vector representing the knots, including the boundary knots.
}
\description{
A helper function that generates \eqn{m + 2} spline knots of the form:
\deqn{
\kappa_{0} = 0 , \kappa_{1} = 1 / (m + 1) , \dots , \kappa_{m} = m / (m + 1) , \kappa_{m + 1} = 1 .
}
The knots are equally spaced with boundary knots \eqn{\kappa_{0} = 0} and \eqn{\kappa_{m + 1} = 1 .}
}
\examples{
generate_knots(3)
}
\references{
Choi, I., Li, B. & Wang, X. (2013). Nonparametric Estimation of Spatial and Space-Time Covariance Function. JABES (Vol. 18, pp. 611-630). 10.1007/s13253-013-0152-z
}
