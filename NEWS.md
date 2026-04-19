# CovEsts 1.0.0

* Initial CRAN submission.

# CovEsts 1.1.0

* Added S3 functionality:
  * All autocovariance estimator functions output S3 objects, called `CovEsts`, and have plot, lines and print methods.
    * `kernel_est` and `splines_est` can output numeric vectors if the argument `estCov` is a numeric vector.
    * `to_pacf` can also output a `CovEsts` object if one is passed.
  * Block bootstrap outputs an S3 object, called `BootEsts`, and has plot, lines and print methods.
  * `to_vario` outputs an S3 object called `VarioEsts` and has plot, lines and print methods.
    * The `VarioEsts` object is created only if a `CovEsts` object is passed as the argument.
  * Transformation functions, such as `make_pd` and the metric functions accept these S3 objects as arguments.
    * `make_pd`, `shrinking` output `CovEsts` objects if they are given as arguments.
* Added parallel computing for `block_bootstrap`, `adjusted_est` and `truncated_est`.
  * This utilises the `parallel` package in provided with `R`.
  
* The following functions was added
  * `normalise_acf`
  
* Several internal functions are no longer user-facing:
  * `H2n`
  * `Xij_mat`
  * `adjusted_spline` 
  * `cyclic_matrix`
  * `generate_knots`
  * `get_taus`
  * `rho_T1`
  * `solve_shrinking` 
  * `solve_spline`
  * `splines_df`

* The following functions were removed.
  * `get_tau`
  * `taper_single`
  * `tapered_single`

* The following functions were deprecated and renamed:
  * `kernel` is now called `kernel_ec`
  * `kernel_symm` is now called `kernel_symm_ec`
  * `window` is now called `window_ec`
  * `window_symm` is now called `window_symm_ec`
