# Most basic call works.
test_that("corrected_est() works", {
  expect_equal(corrected_est(c(1, 2, 3), "gaussian")$acf, c(2/3, 0, -5.398656*10^(-7)))
})

test_that("corrected_est() fails for empty X", {
  expect_error(corrected_est(c(), "gaussian")$acf)
})

test_that("corrected_est() fails for nonnumeric X", {
  expect_error(corrected_est(c(1, 'a', 3), "gaussian")$acf)
  expect_error(corrected_est(c(1, 1i, 3), "gaussian")$acf)
})

test_that("corrected_est() fails for at least one NA in X", {
  expect_error(corrected_est(c(1, NA, 3), "gaussian")$acf)
})

test_that("corrected_est() fails if nonnumeric N_T", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", N_T = 'a')$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", N_T = 1i)$acf)
})

test_that("corrected_est() fails if N_T is NA", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", N_T = NA)$acf)
})

test_that("corrected_est() fails if N_T <= 0", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", N_T = 0)$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", N_T = -0.1)$acf)
})

test_that("corrected_est() fails for nonnumeric meanX", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = 'a')$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = 1i)$acf)
})

test_that("corrected_est() fails for NA meanX", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = NA)$acf)
})

test_that("corrected_est() fails if pd is nonboolean", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", pd = 'TRUE')$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", pd = 1)$acf)
})

test_that("corrected_est() fails if maxLag is nonnumeric", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = 'a')$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = 1i)$acf)
})

test_that("corrected_est() fails if maxLag < 0", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = -1)$acf)
})

test_that("corrected_est() fails if maxLag >= length(X)", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = 3)$acf)
})

test_that("corrected_est() fails for noninteger maxLag", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = 1.5)$acf)
})

test_that("corrected_est() fails if maxLag is NA", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", maxLag = NA)$acf)
})

test_that("corrected_est() fails if meanX's length is not 1", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = c(1, 2))$acf)
})

test_that("corrected_est() fails for nonnumeric meanX", {
    expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = 'a')$acf)
    expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = 1i)$acf)
})

test_that("corrected_est() fails if meanX is NA", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", meanX = NA)$acf)
})

test_that("corrected_est() fails if type is neither 'autocovariance' or 'autocorrelation'", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", type = 'covariance')$acf)
})

test_that("corrected_est() fails if x has at least one NA", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", x = c(1, NA, 3))$acf)
})

test_that("corrected_est() fails for empty x", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", x = c())$acf)
})

test_that("corrected_est() fails for nonnumeric x", {
  expect_error(corrected_est(c(1, 2, 3), "gaussian", x = c(1, 'a', 3))$acf)
  expect_error(corrected_est(c(1, 2, 3), "gaussian", x = c(1, 1i, 3))$acf)
})

test_that("corrected_est() fails for nonboolean custom_kernel", {
  expect_error(corrected_est(c(1, 2, 3), "my_kernel", customKernel = 1)$acf)
  expect_error(corrected_est(c(1, 2, 3), "my_kernel", customKernel = 'TRUE')$acf)
})

# kernel_est
test_that("kernel_est() works", {
  expect_equal(kernel_est(c(1, 2, 3), "gaussian"), c(1, 0.071347986694504816896, 0.000004858790376937838049672957), tolerance = sqrt(.Machine$double.eps))
})

test_that("kernel_est() fails for nonboolean custom_kernel", {
  expect_error(kernel_est(c(1, 2, 3), "my_kernel", customKernel = 1))
  expect_error(kernel_est(c(1, 2, 3), "my_kernel", customKernel = 'TRUE'))
})

test_that("kernel_est() fails for empty estCov", {
expect_error(kernel_est(c(), "gaussian"))
})

test_that("kernel_est() fails for at least one NA in estCov", {
  expect_error(kernel_est(c(1, NA, 3), "gaussian"))
})

test_that("kernel_est() fails for nonnumeric estCov", {
  expect_error(kernel_est(c(1, 'a', 3), "gaussian"))
  expect_error(kernel_est(c(1, 1i, 3), "gaussian"))
})

test_that("kernel_est() fails for nonnumeric N_T", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", N_T = 'a'))
  expect_error(kernel_est(c(1, 2, 3), "gaussian", N_T = 1i))
})

test_that("kernel_est() fails for N_T <= 0", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", N_T = 0))
  expect_error(kernel_est(c(1, 2, 3), "gaussian", N_T = -0.1))
})

test_that("kernel_est() fails for NA N_T", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", N_T = NA))
})

test_that("kernel_est() fails for nonnumeric maxLag", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag='a'))
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag=1i))
})

test_that("kernel_est() fails for maxLag < 0", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag=-1))
})

test_that("kernel_est() fails for maxLag > length(X) - 1", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag=3))
})

test_that("kernel_est() fails for noninteger maxLag", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag=1.5))
})

test_that("kernel_est() fails for NA maxLag", {
  expect_error(kernel_est(c(1, 2, 3), "gaussian", maxLag = NA))
})

# solve_shrinking
test_that("solve_shrinking() fails for nonnumeric par", {
  expect_error(solve_shrinking('a', diag(2), diag(2)))
  expect_error(solve_shrinking(1i, diag(2), diag(2)))
})

test_that("solve_shrinking() fails for NA in par", {
  expect_error(solve_shrinking(NA, diag(2), diag(2)))
})

test_that("solve_shrinking() fails for par not of length 1", {
  expect_error(solve_shrinking(c(0.5, 0.5), diag(2), diag(2)))
})

test_that("solve_shrinking() fails for nonnumeric corr_mat", {
  expect_error(solve_shrinking(0.5, matrix(c(1, 'a', 3, 4), 2), diag(2)))
  expect_error(solve_shrinking(0.5, matrix(c(1, 2i, 3, 4), 2), diag(2)))
})

test_that("solve_shrinking() fails for at least one NA in corr_mat", {
  expect_error(solve_shrinking(0.5, matrix(c(1, NA, 3, 4), 2), diag(2)))
})

test_that("solve_shrinking() fails for nonnumeric target", {
  expect_error(solve_shrinking(0.5, diag(2), matrix(c(1, 'a', 3, 4), 2)))
  expect_error(solve_shrinking(0.5, diag(2), matrix(c(1, 2i, 3, 4), 2)))
})

test_that("solve_shrinking() fails for at least one NA in target", {
  expect_error(solve_shrinking(0.5, diag(2), matrix(c(1, NA, 3, 4), 2)))
})

test_that("solve_shrinking() fails for corr_mat and target having different sizes", {
  expect_error(solve_shrinking(0.5, diag(2), diag(3)))
})

# shrinking
test_that("shrinking() fails for estCov having length less than 1", {
  expect_error(shrinking(c()))
})

test_that("shrinking() fails for nonnumeric estCov", {
  expect_error(shrinking(c(1, 'a', 3)))
  expect_error(shrinking(c(1, 2i, 3)))
})

test_that("shrinking() fails for at least one NA in estCov", {
  expect_error(shrinking(c(1, NA, 3)))
})

test_that("shrinking() fails for nonlogical return_matrix", {
  expect_error(shrinking(c(1, 2, 3), return_matrix = 1))
  expect_error(shrinking(c(1, 2, 3), return_matrix = 'TRUE'))
})

test_that("shrinking() fails for nonmatrix target", {
  expect_error(shrinking(c(1, 2, 3), target = 1))
})

test_that("shrinking() fails for nonnumeric target", {
  expect_error(shrinking(c(1, 2, 3), target = matrix(c(1, 'a', 3, 4), 2)))
  expect_error(shrinking(c(1, 2, 3), target = matrix(c(1, 2i, 3, 4), 2)))
})

test_that("shrinking() fails for at least one NA in target", {
  expect_error(shrinking(c(1, 2, 3), target = matrix(c(1, NA, 3, 4), 2)))
})
