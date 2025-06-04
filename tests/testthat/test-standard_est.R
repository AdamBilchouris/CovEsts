# compute_standard_est works
test_that("compute_standard_est() works", {
  expect_equal(compute_standard_est(c(1, 2, 3)), c(2/3, 0, -1/3))
})

# compute_standard_est fails for empty X
test_that("compute_standard_est() fails for empty X", {
  expect_error(compute_standard_est(c()))
})

# compute_standard_est fails for nonvector X
test_that("compute_standard_est() fails for empty X", {
  expect_error(compute_standard_est(matrix(c(1,2,3,4), 2)))
})

# compute_standard_est fails for nonnumeric X
test_that("compute_standard_est() fails for empty X", {
  expect_error(compute_standard_est(c(1, 'a', 3)))
  expect_error(compute_standard_est(c(1, 1i, 3)))
})

# compute_standard_est fails for nonboolean pd
test_that("compute_standard_est() for nonboolean pd", {
  expect_error(compute_standard_est(c(1, 2, 3), pd='TRUE'))
  expect_error(compute_standard_est(c(1, 2, 3), pd=1))
})

# compute_standard_est fails for nonboolean pd
test_that("compute_standard_est() for nonboolean pd", {
  expect_error(compute_standard_est(c(1, 2, 3), pd='TRUE'))
  expect_error(compute_standard_est(c(1, 2, 3), pd=1))
})

# compute_standard_est fails for maxLag < 0
test_that("compute_standard_est() for maxLag < 0", {
  expect_error(compute_standard_est(c(1, 2, 3), maxLag = -1))
})

# compute_standard_est fails for maxLag > length(X) - 1
test_that("compute_standard_est() for maxLag > length(X) - 1", {
  expect_error(compute_standard_est(c(1, 2, 3), maxLag = 3))
})

# compute_standard_est fails for noninteger maxLag
test_that("compute_standard_est() for noninteger maxLag", {
  expect_error(compute_standard_est(c(1, 2, 3), maxLag = 1.5))
})

# compute_standard_est fails for meanX not of length 1
test_that("compute_standard_est() for meanX not of length 1", {
  expect_error(compute_standard_est(c(1, 2, 3), meanX = c(1, 2)))
})

# compute_standard_est fails for nonnumeric meanX
test_that("compute_standard_est() for nonnumeric meanX", {
  expect_error(compute_standard_est(c(1, 2, 3), meanX = 'a'))
  expect_error(compute_standard_est(c(1, 2, 3), meanX = 1i))
})

# compute_standard_est fails for NA meanX
test_that("compute_standard_est() for NA meanX", {
  expect_error(compute_standard_est(c(1, 2, 3), meanX = NA))
})

# compute_standard_est fails if type is neither autocovariance or autocorrelation
test_that("compute_standard_est() for NA meanX", {
  expect_error(compute_standard_est(c(1, 2, 3), type = 'covariance'))
})
