# standard_est works
test_that("standard_est() works", {
  expect_equal(standard_est(c(1, 2, 3))$acf, c(2/3, 0, -1/3))
})

# standard_est fails for empty X
test_that("standard_est() fails for empty X", {
  expect_error(standard_est(c())$acf)
})

# standard_est fails for nonnumeric X
test_that("standard_est() fails for empty X", {
  expect_error(standard_est(c(1, 'a', 3))$acf)
  expect_error(standard_est(c(1, 1i, 3))$acf)
})

# standard_est fails for nonboolean pd
test_that("standard_est() for nonboolean pd", {
  expect_error(standard_est(c(1, 2, 3), pd='TRUE')$acf)
  expect_error(standard_est(c(1, 2, 3), pd=1)$acf)
})

# standard_est fails for maxLag < 0
test_that("standard_est() for maxLag < 0", {
  expect_error(standard_est(c(1, 2, 3), maxLag = -1)$acf)
})

# standard_est fails for maxLag > length(X) - 1
test_that("standard_est() for maxLag > length(X) - 1", {
  expect_error(standard_est(c(1, 2, 3), maxLag = 3)$acf)
})

# standard_est fails for noninteger maxLag
test_that("standard_est() for noninteger maxLag", {
  expect_error(standard_est(c(1, 2, 3), maxLag = 1.5)$acf)
})

# standard_est fails for meanX not of length 1
test_that("standard_est() for meanX not of length 1", {
  expect_error(standard_est(c(1, 2, 3), meanX = c(1, 2))$acf)
})

# standard_est fails for nonnumeric meanX
test_that("standard_est() for nonnumeric meanX", {
  expect_error(standard_est(c(1, 2, 3), meanX = 'a')$acf)
  expect_error(standard_est(c(1, 2, 3), meanX = 1i)$acf)
})

# standard_est fails for NA meanX
test_that("standard_est() for NA meanX", {
  expect_error(standard_est(c(1, 2, 3), meanX = NA)$acf)
})

# standard_est fails if type is neither autocovariance or autocorrelation
test_that("standard_est() for NA meanX", {
  expect_error(standard_est(c(1, 2, 3), type = 'covariance')$acf)
})

# standard_est fails for empty x
test_that("standard_est() fails for empty x", {
  expect_error(standard_est(c(1, 2, 3), x = c())$acf)
})

# standard_est fails for nonnumeric x
test_that("standard_est() fails nonnumeric x", {
  expect_error(standard_est(c(1, 2, 3), x = c(1, 'a', 3))$acf)
  expect_error(standard_est(c(1, 2, 3), x = c(1, 1i, 3))$acf)
})

# standard_est fails for at least one NA in x
test_that("standard_est() fails for at least one NA in", {
  expect_error(standard_est(c(1, 2, 3), x = c(1, NA, 3))$acf)
})

# to_vario works
test_that("to_vario() works", {
  expect_equal(to_vario(c(1, 0.5, 0)), c(0, 0.5, 1))
})

# to_vario fails if estCov is empty
test_that("to_vario() fails if estCov is empty", {
  expect_error(to_vario(c()))
})

# to_vario fails if at least one value in estCov is nonnumeric
test_that("to_vario() fails if at least one value in estCov is nonnumeric", {
  expect_error(to_vario(c(1, 'a', 3)))
  expect_error(to_vario(c(1, 1i, 3)))
})

# to_vario fails if estCov has at least one NA
test_that("to_vario() fails f estCov has at least one NA", {
  expect_error(to_vario(c(1, NA, 3)))
})

# to_pacf works
test_that("to_pacf() works", {
  expect_equal(to_pacf(c(1, 0.5, 0)), c(0.5, -1/3))
})

# to_pacf fails if estCov is empty
test_that("to_pacf() fails if estCov is empty", {
  expect_error(to_pacf(c()))
})

# to_pacf fails if at least one value in estCov is nonnumeric
test_that("to_pacf() fails if at least one value in estCov is nonnumeric", {
  expect_error(to_pacf(c(1, 'a', 3))$acf)
  expect_error(to_pacf(c(1, 1i, 3))$acf)
})

# to_pacf fails if estCov has at least one NA
test_that("to_pacf() fails if estCov has at least one NA", {
  expect_error(to_pacf(c(1, NA, 3))$vario)
})
