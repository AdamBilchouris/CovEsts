# CovEsts: Nonparametric Autocovariance Estimation and Analysis <img src="man/figures/logo.png" align="right" height="120" alt="" />

<!-- badges: start -->
<!-- badges: end -->


The estimation of the autocovariance function is an important problem for both theory and applied fields.


Several nonparametric autocovariance function estimators are implemented in this package, each with their own statistical properties.
Having multiple estimators allows one to obtain estimates with different behaviour, which can lead to more appropriate estimates.
The implemented estimators are for the one-dimensional stochastic processes.


## Publication
Bilchouris, A., & Olenko, A. (2025). On Nonparametric Estimation of Covariogram. Austrian Journal of Statistics, 54(1), 112–137. https://doi.org/10.17713/ajs.v54i1.1975


## Features
`CovEsts` provides eight different estimators of the autocovariance function, with each providing theoretical and practical benefits.
Briefly, the estimators are

* the classic estimators,
* two kernel regression estimators,
* an estimator adjusting for the edge effect,
* an estimator utilising splines,
* kernel correction of estimators,
* computing the semivariogram from an autocovariance function estimate, and
* computing the partial autocorelation function from an autocovariance function estimate.

The package also provides several general functions, which are
  
* the forward and inverse one-dimensional discrete cosine transforms,
* making a function positive-definite,
* linear shrinking, and
* computing the nearest positive-definite matrix.

## Installation

You can install the development version of CovEsts from [GitHub](https://github.com/AdamBilchouris/CovEsts) with:

``` r
# install.packages("devtools")
devtools::install_github("AdamBilchouris/CovEsts")
```

## Example

Below is a basic use case of the package. 

``` r
library(CovEsts)
## basic example code
X <- rnorm(100)
standard_est(X, maxLag = 20)
plot(standard_est(X, maxLag = 20))
```

