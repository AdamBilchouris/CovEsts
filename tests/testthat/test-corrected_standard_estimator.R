# Most basic call works.
test_that("corrected_standard_estimator() works for upperTau=2, N=3, positive-definite, nonzero mean.", {
  expect_equal(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian"), c(2/3, 0, -5.398656*10^(-7)))
})

# Fails for X not being a vector.
test_that("corrected_standard_estimator() fails if X is not a vector", {
  expect_error(corrected_standard_estimator(1, 2, "gaussian"))
  expect_error(corrected_standard_estimator(matrix(rnorm(4), 2, 2), 2, "gaussian"))
})

# An unknown kernel should error.
test_that("corrected_standard_estimator() fails for an unknown kernel", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "test_kernel"))
})

# Length should be greater than zero.
test_that("corrected_standard_estimator() fails for an empty process", {
  expect_error(corrected_standard_estimator(c(), 2, "gaussian"))
})

# X should be a numeric vector,
test_that("corrected_standard_estimator() fails for a nonnumeric X", {
  expect_error(corrected_standard_estimator(c(1i, 'a'), 2, "gaussian"))
})

# upperTau must be at least 0
test_that("corrected_standard_estimator() fails for a negative upperTau", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), -1, "gaussian"))
})

# upperTau must be less than the length of the process.
test_that("corrected_standard_estimator() fails for an upperTau greater than the length of the process", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 3, "gaussian"))
})

# pd (positive-definite) must be a boolean.
test_that("corrected_standard_estimator() fails for a nonboolean pd", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", pd=1))
})

# type can only be covariance or correlation.
test_that("corrected_standard_estimator() fails for type being neither covariance nor correlation", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", type="neither"))
})

# customKernel must be a boolean.
test_that("corrected_standard_estimator() fails for custom kernel being a nonboolean", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", customKernel=1))
})

# N is different from the length of the process.
test_that("corrected_standard_estimator() fails N being different from the length of the process", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", N=1))
})

# N_T fails for a nonnumeric case and being not being greater than zero.
test_that("corrected_standard_estimator() fails for N_T being nonnumeric and not being greater than zero", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", N_T='a'))
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", N=0))
})

# Fails for nonnumeric meanX
test_that("corrected_standard_estimator() fails for meanX being nonnumeric", {
  expect_error(corrected_standard_estimator(c(1, 2, 3), 2, "gaussian", meanX=1i))
})

# Works for a custom kernel.
# Figure out how to get this to work.
# test_that("corrected_standard_estimator() works for a custom kernel", {
#   my_kernel <- function(x, theta, params) {
#     stopifnot(theta > 0, length(x) >= 1, all(x >= 0))
#     return(sapply(x, function(t) ifelse(t == 0, 1, ifelse(t == Inf, 0,(sin((t^params[1]) / theta) / ((t^params[1]) / theta)) * cos((t^params[2]) / theta)))))
#   }
#   expect_equal(corrected_standard_estimator(c(1, 2, 3), 2, "my_kernel", kernel_params=c(2, 0.25), customKernel = TRUE), c(0.66666667, 0, 0.01180484))
# })
