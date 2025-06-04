# Most basic call works.
test_that("compute_corrected_standard_est() works", {
  expect_equal(compute_corrected_standard_est(c(1, 2, 3), "gaussian"), c(2/3, 0, -5.398656*10^(-7)))
})

test_that("compute_corrected_standard_est() fails for empty X", {
  expect_error(compute_corrected_standard_est(c(), "gaussian"))
})

test_that("compute_corrected_standard_est() fails if X is not a vector", {
  expect_error(compute_corrected_standard_est(matrix(c(1, 2, 3, 4), 2, 2), "gaussian"))
})

test_that("compute_corrected_standard_est() fails for nonnumeric X", {
  expect_error(compute_corrected_standard_est(c(1, 'a', 3), "gaussian"))
  expect_error(compute_corrected_standard_est(c(1, 1i, 3), "gaussian"))
})

test_that("compute_corrected_standard_est() fails if X is not a vector", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", N_T = 'a'))
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", N_T = 1i))
})

test_that("compute_corrected_standard_est() fails if N_T <= 0", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", N_T = 0))
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", N_T = -0.1))
})

test_that("compute_corrected_standard_est() fails for nonnumeric meanX", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", meanX = 'a'))
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", meanX = 1i))
})

test_that("compute_corrected_standard_est() fails if pd is nonboolean", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", pd = 'TRUE'))
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", pd = 1))
})

test_that("compute_corrected_standard_est() fails if maxLag is nonnumeric", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = 'a'))
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = 1i))
})

test_that("compute_corrected_standard_est() fails if maxLag < 0", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = -1))
})

test_that("compute_corrected_standard_est() fails if maxLag >= length(X)", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = 3))
})

test_that("compute_corrected_standard_est() fails for noninteger maxLag", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = 1.5))
})

test_that("compute_corrected_standard_est() fails if meanX's length is not 1", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", meanX = c(1, 2)))
})

test_that("compute_corrected_standard_est() fails for nonnumeric meanX", {
    expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", meanX = 'a'))
    expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", meanX = 1i))
})

test_that("compute_corrected_standard_est() fails if meanX is NA", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", maxLag = NA))
})

test_that("compute_corrected_standard_est() fails if type is neither 'autocovariance' or 'autocorrelation'", {
  expect_error(compute_corrected_standard_est(c(1, 2, 3), "gaussian", type = 'covariance'))
})

# compute_kernel_corrected_est
test_that("compute_kernel_corrected_est() works", {
  expect_equal(compute_kernel_corrected_est(c(1, 2, 3), "gaussian"), c(1, 0.071347986694504816896, 0.000004858790376937838049672957), tolerance = sqrt(.Machine$double.eps))
})

test_that("compute_kernel_corrected_est() fails for nonboolean custom_kernel", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "my_kernel", customKernel = 1))
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "my_kernel", customKernel = 'TRUE'))
})

test_that("compute_kernel_corrected_est() fails for empty cov", {
expect_error(compute_kernel_corrected_est(c(), "gaussian"))
})

test_that("compute_kernel_corrected_est() fails for nonvector cov", {
  expect_error(compute_kernel_corrected_est(matrix(c(1, 2, 3, 4), 2), "gaussian"))
})

test_that("compute_kernel_corrected_est() fails for nonnumeric cov", {
  expect_error(compute_kernel_corrected_est(c(1, 'a', 3), "gaussian"))
  expect_error(compute_kernel_corrected_est(c(1, 1i, 3), "gaussian"))
})

test_that("compute_kernel_corrected_est() fails for nonnumeric N_T", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", N_T = 'a'))
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", N_T = 1i))
})

test_that("compute_kernel_corrected_est() fails for N_T <= 0", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", N_T = 0))
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", N_T = -0.1))
})

test_that("compute_kernel_corrected_est() fails for nonnumeric maxLag", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", maxLag='a'))
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", maxLag=1i))
})

test_that("compute_kernel_corrected_est() fails for maxLag < 0", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", maxLag=-1))
})

test_that("compute_kernel_corrected_est() fails for maxLag > length(X) - 1", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", maxLag=3))
})

test_that("compute_kernel_corrected_est() fails for noninteger maxLag", {
  expect_error(compute_kernel_corrected_est(c(1, 2, 3), "gaussian", maxLag=1.5))
})

