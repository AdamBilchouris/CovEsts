
# CovEsts: Nonparametric Autocovariance Estimation and Estimator Adjustment

<!-- badges: start -->
<!-- badges: end -->

This package provides estimates of the autocovariance function for one-dimensional processes.

# Publication
Bilchouris, A., & Olenko, A. (2025). On Nonparametric Estimation of Covariogram. Austrian Journal of Statistics, 54(1), 112–137. https://doi.org/10.17713/ajs.v54i1.1975

## Installation

You can install the development version of CovEst from [GitHub](https://github.com/AdamBilchouris/CovEsts) with:

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
compute_standard_est(X, 20)
plot(compute_standard_est(X, 20))
```

