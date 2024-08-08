# Positive-definite
test_that("standard_est_single() works for tau = 1, N = 3, positive-definite, nonzero mean.", {
  expect_equal(standard_est_single(c(1, 2, 3), 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE), 0)
})

test_that("standard_est_single() fails for empty X, tau=0, N = 3, positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(), 0, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE))
})

test_that("standard_est_single() fails for tau = -1, N = 3, positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), -1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE))
})

test_that("standard_est_single() fails for tau >= N, N = 3, positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), length(c(1, 2, 3)), N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = TRUE))
})

test_that("standard_est_single() fails for N > length(X), tau=0, N = 4, positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), 0, N = 4, meanX = mean(c(1, 2, 3)), pd = TRUE))
})

test_that("standard_est_single() fails for pd is neither TRUE or FALSE, tau=0, N = 4, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), 0, N = length(c(1, 2 ,3)), meanX = mean(c(1, 2, 3)), pd = 'TRUE'))
})

# Nonpositive-definite
test_that("standard_est_single() works for tau = 1, N = 3, not positive-definite, nonzero mean.", {
  expect_equal(standard_est_single(c(1, 2, 3), 1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE), 0)
})

test_that("standard_est_single() fails for empty X, tau=0, N = 3, not positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(), 0, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE))
})

test_that("standard_est_single() fails for tau = -1, N = 3, not positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), -1, N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE))
})

test_that("standard_est_single() fails for tau >= N, N = 3, not positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), length(c(1, 2, 3)), N = length(c(1, 2, 3)), meanX = mean(c(1, 2, 3)), pd = FALSE))
})

test_that("standard_est_single() fails for N > length(X), tau=0, N = 4, not positive-definite, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), 0, N = 4, meanX = mean(c(1, 2, 3)), pd = FALSE))
})

test_that("standard_est_single() fails for pd is neither TRUE or FALSE, tau=0, N = 4, nonzero mean.", {
  expect_error(standard_est_single(c(1, 2, 3), 0, N = length(c(1, 2 ,3)), meanX = mean(c(1, 2, 3)), pd = 'FALSE'))
})
