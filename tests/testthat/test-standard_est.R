# Positive-definite
test_that("standard_est() works for upperTau = N - 1, N = 3, positive-definite, nonzero mean, covariance.", {
  expect_equal(standard_est(c(1, 2, 3), length(c(1, 2, 3)) - 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='covariance'), c(2/3, 0, -1/3))
})

test_that("standard_est() works for upperTau = N - 1, N = 3, positive-definite, nonzero mean, correlation", {
  expect_equal(standard_est(c(1, 2, 3), length(c(1, 2, 3)) - 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='correlation'), c(1, 0, -1/2))
})

test_that("standard_est() fails for empty X, upperTau=2, N = 3, positive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(), 2, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='covariance'))
})

test_that("standard_est() fails for upperTau = -1, N = 3, positive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), -1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='covariance'))
})

test_that("standard_est() fails for upperTau >= N, N = 3, positive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), 3, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='covariance'))
})

test_that("standard_est() fails for N > length(X), upperTau=2, N = 4, positive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), 2, N = 4, meanX = mean(c(1, 2, 3)), pd = TRUE, type='covariance'))
})

test_that("standard_est() fails for pd is neither TRUE or FALSE, tau=0, N = 4, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), 0, N = length(c(1, 2 ,3)), meanX = mean(c(1, 2, 3)), pd = 'TRUE', type='covariance'))
})

test_that("standard_est() fails for type is neither 'covariance' or 'correlation', upperTau=2, N = 4, nonzero mean.", {
  expect_error(standard_est(c(1, 2, 3), 2, N = length(c(1, 2 ,3)), meanX = mean(c(1, 2, 3)), pd = TRUE, type='cov'))
})

# Nonpositive-definite
test_that("standard_est() works for upperTau = N - 1, N = 3, nonpositive-definite, nonzero mean, covariance.", {
  expect_equal(standard_est(c(1, 2, 3), length(c(1, 2, 3)) - 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='covariance'), c(2/3, 0, -1))
})

test_that("standard_est() works for upperTau = N - 1, N = 3, nonpositive-definite, nonzero mean, correlation", {
  expect_equal(standard_est(c(1, 2, 3), length(c(1, 2, 3)) - 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='correlation'), c(1, 0, -3/2))
})

test_that("standard_est() fails for empty X, upperTau=2, N = 3, nonpositive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(), 2, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='covariance'))
})

test_that("standard_est() fails for upperTau = -1, N = 3, nonpositive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), -1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='covariance'))
})

test_that("standard_est() fails for upperTau >= N, N = 3, nonpositive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), 3, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='covariance'))
})

test_that("standard_est() fails for N > length(X), upperTau=2, N = 4, nonpositive-definite, nonzero mean, covariance.", {
  expect_error(standard_est(c(1, 2, 3), 2, N = 4, meanX = mean(c(1, 2, 3)), pd = FALSE, type='covariance'))
})

test_that("standard_est() fails for type is neither 'covariance' or 'correlation', upperTau=2, N = 4, nonzero mean.", {
  expect_error(standard_est(c(1, 2, 3), 2, N = length(c(1, 2 ,3)), meanX = mean(c(1, 2, 3)), pd = FALSE, type='cov'))
})
