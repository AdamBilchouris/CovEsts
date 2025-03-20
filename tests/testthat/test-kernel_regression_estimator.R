# Xij_mat
test_that("Xij_mat works", {
  expect_equal(Xij_mat(c(1, 2)), matrix(c(0.25, -0.25, -0.25, 0.25), ncol=2, byrow=T))
})

test_that("Xij_mat fails for an empty process", {
  expect_error(Xij_mat(c()))
})

test_that("Xij_mat fails for nonnumeric process", {
  expect_error(Xij_mat(c('a')))
  expect_error(Xij_mat(c(1i)))
})

test_that("Xij_mat fails for a process with at least one NA", {
  expect_error(Xij_mat(c(1, NA, 2)))
})

# dct_1d
test_that("dct_1d works", {
  expect_equal(dct_1d(c(1, 2, 3)), c(6.000000, -1.732051, 0.000000), tolerance = 1e-6)
})

test_that("dct_1d fails for nonnumeric X", {
  expect_error(dct_1d(c('a')))
  expect_error(dct_1d(c(1i)))
})

test_that("dct_1d fails for empty X", {
  expect_error(dct_1d(c()))
})

test_that("dct_1d fails for at least one NA", {
  expect_error(dct_1d(c(1, NA, 3)))
})

# idct_1d
test_that("idct_1d works", {
  expect_equal(idct_1d(c(6.000000, -1.732051, 0.00000)), c(1, 2, 3), tolerance = 1e-6)
})

test_that("idct_1d fails for nonnumeric X", {
  expect_error(idct_1d(c('a')))
  expect_error(idct_1d(c(1i)))
})

test_that("idct_1d fails for empty X", {
  expect_error(idct_1d(c()))
})

test_that("idct_1d fails for at least one NA", {
  expect_error(idct_1d(c(1, NA, 3)))
})

# rho_T1
test_that("rho_T1 works", {
  expect_equal(rho_T1(test_x, test_meanX, 2, 0.01, test_xij), -1)
})

test_that("rho_T1 faisl for nonnumeric x", {
  expect_error(rho_T1(c('a'), test_meanX, 2, 0.01, test_xij))
  expect_error(rho_T1(c(1i), test_meanX, 2, 0.01, test_xij))
})

test_that("rho_T1 fails for empty x", {
  expect_error(rho_T1(c(), test_meanX, 2, 0.01, test_xij))
})

test_that("rho_T1 fails for at least one NA in x", {
  expect_error(rho_T1(c(1, NA, 2), test_meanX, 2, 0.01, test_xij))
})

test_that("rho_T1 fails for meanX of length not equal to 1", {
  expect_error(rho_T1(test_x, c(1, 2), 2, 0.01, test_xij))
  expect_error(rho_T1(test_x, c(), 2, 0.01, test_xij))
})

test_that("rho_T1 fails for nonnumeric meanX", {
  expect_error(rho_T1(test_x, 'a', 2, 0.01, test_xij))
  expect_error(rho_T1(test_x, 1i, 2, 0.01, test_xij))
})

test_that("rho_T1 fails for NA meanX", {
  expect_error(rho_T1(test_x, NA, 2, 0.01, test_xij))
})

test_that("rho_T1 fails for T1 of length not equal to 1", {
  expect_error(rho_T1(test_x, test_meanX, c(1, 2), 0.01, test_xij))
  expect_error(rho_T1(test_x, test_meanX, c(), 0.01, test_xij))
})

test_that("rho_T1 fails for nonnumeric T1", {
  expect_error(rho_T1(test_x, test_meanX, 'a', 0.01, test_xij))
  expect_error(rho_T1(test_x, test_meanX, 1i, 0.01, test_xij))
})

test_that("rho_T1 fails for NA T1", {
  expect_error(rho_T1(test_x, test_meanX, NA, 0.01, test_xij))
})

test_that("rho_T1 fails for nonpositive T1", {
  expect_error(rho_T1(test_x, test_meanX, 0, 0.01, test_xij))
  expect_error(rho_T1(test_x, test_meanX, -1, 0.01, test_xij))
})

test_that("rho_T1 fails for h of length not equal to 1", {
  expect_error(rho_T1(test_x, test_meanX, 2, c(1, 2), test_xij))
  expect_error(rho_T1(test_x, test_meanX, 2, c(), test_xij))
})

test_that("rho_T1 fails for nonnumeric h", {
  expect_error(rho_T1(test_x, test_meanX, 2, 'a', test_xij))
  expect_error(rho_T1(test_x, test_meanX, 2, 1i, test_xij))
})

test_that("rho_T1 fails for NA h", {
  expect_error(rho_T1(test_x, test_meanX, 2, NA, test_xij))
})

test_that("rho_T1 fails for nonpositive h", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij))
  expect_error(rho_T1(test_x, test_meanX, 2, -1, test_xij))
})

test_that("rho_T1 fails for nonnumeric Xij_mat", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0.01, matrix(c('a', 1, 2, 1), ncol=2)))
  expect_error(rho_T1(test_x, test_meanX, 2, 0.01, matrix(c(1i, 1, 2, 1), ncol=2)))
})

test_that("rho_T1 fails for nonmatrix Xij_mat", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0.01, 1))
  expect_error(rho_T1(test_x, test_meanX, 2, 0.01, c(1, 2)))
  expect_error(rho_T1(test_x, test_meanX, 2, 0.01, list(a=1, b=2)))
})

test_that("rho_T1 fails for any NA in Xij_mat", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0, matrix(c(1, 1, NA, 1), ncol=2)))
})

test_that("rho_T1 fails for any nonboolean custom_kernel", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij, custom_kernel = 1))
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij, custom_kernel = 'TRUE'))
})

test_that("rho_T1 fails for any nonboolean custom_kernel", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij, custom_kernel = 1))
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij, custom_kernel = 'TRUE'))
})

test_that("rho_T1 fails for kernel_name not in 'gaussian' 'wave', 'rational_quadratic', 'bessel_j'", {
  expect_error(rho_T1(test_x, test_meanX, 2, 0, test_xij, kernel_name = "123", custom_kernel = TRUE))
})

# compute_truncated_est
test_that("compute_truncated_est works", {
  expect_equal(compute_truncated_est(test_X, test_x, c(1, 2, 3), 2, 3, 0.01, "gaussian", meanX = test_meanX), c(0, -1, 0))
})

test_that("compute_truncated_est fails for nonnumeric X", {
  expect_error(compute_truncated_est(c(1, 'a', 3), test_x, 2, 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(c(1, 1i, 3), test_x, 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for empty X", {
  expect_error(compute_truncated_est(c(), test_x, 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for at least one NA in X", {
  expect_error(compute_truncated_est(c(1, NA, 3), test_x, 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for empty x", {
  expect_error(compute_truncated_est(test_X, c(), 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric x", {
  expect_error(compute_truncated_est(test_X, c(1, 'a', 3), 2, 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, c(1, 1i, 3), 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for at least one NA in x", {
  expect_error(compute_truncated_est(test_X, c(1, NA, 3), 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for meanX of length not equal to 1", {
  expect_error(compute_truncated_est(test_X, test_x, c(1, 2), 2, 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, c(), 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric meanX", {
  expect_error(compute_truncated_est(test_X, test_x, 'a', 2, 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 1i, 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for NA meanX", {
  expect_error(compute_truncated_est(test_X, test_x, NA, 2, 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for empty t", {
  expect_error(compute_truncated_est(test_X, test_x, c(), 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric t", {
  expect_error(compute_truncated_est(test_X, test_x, 'a', 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 1i, 2, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, c(1, 'a', 3), 2, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for T1 of length not equal to 1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, c(1, 2), 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, c(), 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric T1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 'a', 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 1i, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for NA T1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, NA, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for T1 <= 0", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 0, 3, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, -0.5, 3, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for T2 of length not equal to 1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, c(1, 2), 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, c(), 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric T2", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 'a', 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 1i, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for NA T2", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, NA, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for T2 <= T1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 1, 0.01, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 2, 0.01, "gaussian"))
})

test_that("compute_truncated_est fails for h of length not equal to 1", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, c(1, 2), "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, c(), "gaussian"))
})

test_that("compute_truncated_est fails for nonnumeric h", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 'a', "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 1i, "gaussian"))
})

test_that("compute_truncated_est fails for nonpositive h", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 0, "gaussian"))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, -0.01, "gaussian"))
})

test_that("compute_truncated_est fails for nonboolean custom_kernel", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 0.01, "gaussian", custom_kernel = 1))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 0.01, "gaussian", custom_kernel = 'TRUE'))
})

test_that("compute_truncated_est fails for nonboolean pd", {
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 0.01, "gaussian", pd = 1))
  expect_error(compute_truncated_est(test_X, test_x, 2, 2, 3, 0.01, "gaussian", pd = 'TRUE'))
})

# compute_adjusted_est_est
test_that("compute_adjusted_est works", {
  expect_equal(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, "gaussian", meanX = test_meanX), c(2/3, -1/3, -1/3))
})

test_that("compute_adjusted_est fails for nonnumeric X", {
  expect_error(compute_adjusted_est(c(1, 'a', 3), test_x, c(1, 2, 3), 0.01, "gaussian"))
  expect_error(compute_adjusted_est(c(1, 1i, 3), test_x, c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for empty X", {
  expect_error(compute_adjusted_est(c(), test_x, c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for at least one NA in X", {
  expect_error(compute_adjusted_est(c(1, NA, 3), test_x, c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for nonnumeric x", {
  expect_error(compute_adjusted_est(test_X, c(1, 'a', 3), c(1, 2, 3), 0.01, "gaussian"))
  expect_error(compute_adjusted_est(test_X, c(1, 1i, 3), c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for empty x", {
  expect_error(compute_adjusted_est(test_X, c(), c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for at least one NA in x", {
  expect_error(compute_adjusted_est(test_X, c(1, NA, 3), c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for meanX of length not equal to 1", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2), c(1, 2, 3), 0.01, "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, c(), c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for nonnumeric meanX", {
  expect_error(compute_adjusted_est(test_X, test_x, 'a', c(1, 2, 3), 0.01, "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, 1i, c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for NA meanX", {
  expect_error(compute_adjusted_est(test_X, test_x, NA, c(1, 2, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for at least one NA in t", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, NA, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for nonnumeric t", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 'a', 3), 0.01, "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 1i, 3), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for empty t", {
  expect_error(compute_adjusted_est(test_X, test_x, c(), 0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for h of length not equal to 1", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), c(1, 2), "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), c(), "gaussian"))
})

test_that("compute_adjusted_est fails for nonnumeric h", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 'a', "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 1i, "gaussian"))
})

test_that("compute_adjusted_est fails for h <= 0", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0, "gaussian"))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), -0.01, "gaussian"))
})

test_that("compute_adjusted_est fails for nonboolean custom_kernel", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, "gaussian", custom_kernel = 1))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, "gaussian", custom_kernel = 'TRUE'))
})

test_that("compute_adjusted_est fails for nonboolean pd", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, "gaussian", pd = 1))
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, "gaussian", pd = 'TRUE'))
})

test_that("compute_adjusted_est fails for kernel_name not in 'gaussian' 'wave', 'rational_quadratic', 'bessel_j'", {
  expect_error(compute_adjusted_est(test_X, test_x, c(1, 2, 3), 0.01, kernel_name = "123"))
})

# make_pd
test_that("make_pd works", {
  expect_equal(make_pd(c(1, 0, 3), TRUE), c(2, 0, 2))
  expect_equal(make_pd(c(1, 0, 3), FALSE), c(4/3, 4/3, 4/3))
})


test_that("make_pd fails for nonnumeric x", {
  expect_error(make_pd(c(1, 'a', 3), TRUE))
  expect_error(make_pd(c(1, 1i, 3), TRUE))
})

test_that("make_pd fails for empty x", {
  expect_error(make_pd(c(), TRUE))
})

test_that("make_pd fails for at least one NA in x", {
  expect_error(make_pd(c(1, NA, 3), TRUE))
})

test_that("make_pd fails for nonboolean method.1", {
  expect_error(make_pd(c(1, 2, 3), 1))
  expect_error(make_pd(c(1, 2, 3), 'TRUE'))
})
