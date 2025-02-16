
# CovEsts

<!-- badges: start -->
<!-- badges: end -->

This package provides estimates of the autocovariance function for one-dimensional processes.

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

