# generate_knots()
test_that("generate_knots() works", {
  expect_equal(generate_knots(3), c(0, 1/4, 1/2, 3/4, 1))
})

test_that("generate_knots() fails for nonnumeric m", {
  expect_error(generate_knots('a'))
  expect_error(generate_knots(1i))
})

test_that("generate_knots() fails for m <= 0", {
  expect_error(generate_knots(0))
})

test_that("generate_knots() fails for noninteger m", {
  expect_error(generate_knots(1.1))
})

# get_tau
test_that("get_tau() works", {
  expect_equal(get_tau(1, test_p, test_m, test_kVec), 1/4)
})

test_that("get_tau() fails for nonnumeric p", {
  expect_error(get_tau(1, 'a', test_m, test_kVec))
  expect_error(get_tau(1, 1i, test_m, test_kVec))
})

test_that("get_tau() fails for p < 0", {
  expect_error(get_tau(1, -0.1, test_m, test_kVec))
})

test_that("get_tau() fails for noninteger p", {
  expect_error(get_tau(1, 0.1, test_m, test_kVec))
})

test_that("get_tau() fails for nonnumeric m", {
  expect_error(get_tau(1, test_p, 'a', test_kVec))
  expect_error(get_tau(1, test_p, 1i, test_kVec))
})

test_that("get_tau() fails for m < 0", {
  expect_error(get_tau(1, test_p, -0.1, test_kVec))
})

test_that("get_tau() fails for noninteger m", {
  expect_error(get_tau(1, test_p, 0.1, test_kVec))
})

test_that("get_tau() fails for nonvector kVec", {
  expect_error(get_tau(1, test_p, test_m, matrix(c(1, 2, 3, 4), nrow=2)))
})

test_that("get_tau() fails for empty kVec", {
  expect_error(get_tau(1, test_p, test_m, c()))
})

test_that("get_tau() fails for any NA in kVec", {
  expect_error(get_tau(1, test_p, test_m, c(0, NA, 1)))
})

test_that("get_tau() fails if any elements of kVec are not between 0 and 1", {
  expect_error(get_tau(1, test_p, test_m, c(0, 0.5, 1.1)))
})

test_that("get_tau() returns NA if i > -p or i > m + p + 1)", {
  expect_equal(get_tau(-test_p - 1, test_p, test_m, test_kVec), NA)
  expect_equal(get_tau(test_p + test_m + 2, test_p, test_m, test_kVec), NA)
})

test_that("get_all_tau() works", {
  expect_equal(get_all_tau(test_p, test_m), list('-2' = -1/2, '-1' = -1/4, '0' = 0, '1' = 1/4, '2'= 1/2, '3' = 3/4, '4' = 1, '5' = 5/4, '6' = 3/2))
})

test_that("get_all_tau() fails for nonnumeric p", {
  expect_error(get_all_tau('a', test_m))
  expect_error(get_all_tau(1i, test_m))
})

test_that("get_all_tau() fails for p < 0", {
  expect_error(get_all_tau(-0.1, test_m))
})

test_that("get_all_tau() fails for noninteger p", {
  expect_error(get_all_tau(1.1, test_m))
})

test_that("get_all_tau() fails for nonnumeric m", {
  expect_error(get_all_tau(test_p, 'a'))
  expect_error(get_all_tau(test_p, 1i))
})

test_that("get_all_tau() fails for m < 0", {
  expect_error(get_all_tau(test_p, -0.1))
})

test_that("get_all_tau() fails for noninteger m", {
  expect_error(get_all_tau(test_p, 1.1))
})

# f_j_l
test_that("f_j_l() works", {
  expect_equal(f_j_l(1, 1, 1, test_p, test_m, test_taus), 1/24)
})

test_that("f_j_l() fails for nonnumeric x", {
  expect_error(f_j_l('a', 1, 1, test_p, test_m, test_taus))
  expect_error(f_j_l(1i, 1, 1, test_p, test_m, test_taus))
})

test_that("f_j_l() fails for noninteger j", {
  expect_error(f_j_l(1, 1.1, 1, test_p, test_m, test_taus))
})

test_that("f_j_l() fails for j <= 0", {
  expect_error(f_j_l(1, -1, 1, test_p, test_m, test_taus))
})

test_that("f_j_l() fails for noninteger l", {
  expect_error(f_j_l(1, 1, 1.1, test_p, test_m, test_taus))
})

test_that("f_j_l() fails for nonnumeric l", {
  expect_error(f_j_l(1, 1, 'a', test_p, test_m, test_taus))
  expect_error(f_j_l(1, 1, 1i, test_p, test_m, test_taus))
})

test_that("f_j_l() fails for nonnumeric p", {
  expect_error(f_j_l(1, 1, 1, 'a', test_m, test_taus))
  expect_error(f_j_l(1, 1, 1, 1i, test_m, test_taus))
})

test_that("f_j_l() fails for p < 0", {
  expect_error(f_j_l(1, 1, 1, -1, test_m, test_taus))
})

test_that("f_j_l() fails for noninteger p", {
  expect_error(f_j_l(1, 1, 1, 1.1, test_m, test_taus))
})

test_that("f_j_l() fails for nonnumeric m", {
  expect_error(f_j_l(1, 1, 1, test_p, 'a', test_taus))
  expect_error(f_j_l(1, 1, 1, test_p, 1i, test_taus))
})

test_that("f_j_l() fails for m < 0", {
  expect_error(f_j_l(1, 1, 1, test_p, -1, test_taus))
})

test_that("f_j_l() fails for noninteger m", {
  expect_error(f_j_l(1, 1, 1, test_p, 1.1, test_taus))
})

test_that("f_j_l() fails for nonvector taus", {
  expect_error(f_j_l(1, 1, 1, test_p, test_m, matrix(c(1, 2, 3, 4), 2)))
})

test_that("f_j_l() fails for empty taus", {
  expect_error(f_j_l(1, 1, 1, test_p, test_m, c()))
})

test_that("f_j_l() fails for at least one NA in taus", {
  expect_error(f_j_l(1, 1, 1, test_p, test_m, c(1, NA, 2)))
})

# get_splines_df
test_that("get_splines_df() fails for nonnumeric x", {
  expect_equal(get_splines_df(1, 1, 1, list('-1' = -0.5, '0' = 0, '1' = 0.5, '2' = 1, '3' = 1.5)),
               data.frame(lags = 1, j1 = 1/4, j2 = 3/4))
})
test_that("get_splines_df() fails for nonnumeric x", {
  expect_error(get_splines_df(c(1, 'a', 3), test_p, test_m, test_taus))
  expect_error(get_splines_df(c(1, 1i, 3), test_p, test_m, test_taus))
})

test_that("get_splines_df() fails for empty x", {
  expect_error(get_splines_df(c(), test_p, test_m, test_taus))
})

test_that("get_splines_df() fails for at least one NA in x", {
  expect_error(get_splines_df(c(1, NA, 3), test_p, test_m, test_taus))
})

test_that("get_splines_df() fails for nonnumeric p", {
  expect_error(get_splines_df(test_x, 'a', test_m, test_taus))
  expect_error(get_splines_df(test_x, 1i, test_m, test_taus))
})

test_that("get_splines_df() fails for p < 0", {
  expect_error(get_splines_df(test_x, -1, test_m, test_taus))
})

test_that("get_splines_df() fails for nonninteger p", {
  expect_error(get_splines_df(test_x, 1.1, test_m, test_taus))
})

test_that("get_splines_df() fails for nonnumeric m", {
  expect_error(get_splines_df(test_x, test_p, 'a', test_m, test_taus))
  expect_error(get_splines_df(test_x, test_p, 1i, test_taus))
})

test_that("get_splines_df() fails for m < 0", {
  expect_error(get_splines_df(test_x, test_p, -1, test_taus))
})

test_that("get_splines_df() fails for nonninteger m", {
  expect_error(get_splines_df(test_x, test_p, 1.1, test_taus))
})

test_that("get_splines_df() fails for nonvector taus", {
  expect_error(get_splines_df(test_x, test_p, test_m, matrix(c(1, 2, 3, 4), 2)))
})

test_that("get_splines_df() fails for empty taus", {
  expect_error(get_splines_df(test_x, test_p, test_m, c()))
})

test_that("get_splines_df() fails for at least one NA in taus", {
  expect_error(get_splines_df(test_x, test_p, test_m, c(1, NA, 2)))
})

# compute_splines_est
#test_that("compute_splines_est() works", {
#  expect_equal(compute_splines_est(test_X, test_x, 2, test_estCov, test_p, test_m), c(0.666666184, 0.002083332), tolerance = sqrt(.Machine$double.eps))
#})
# hard to test this function as it is random.

test_that("compute_splines_est() fails for nonnumeric X", {
  expect_error(compute_splines_est(c(1, 'a', 3), test_x, test_estCov, test_p, test_m, maxLag = 2))
  expect_error(compute_splines_est(c(1, 1i, 3), test_x, test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonvector X", {
  expect_error(compute_splines_est(matrix(c(1, 2, 3, 4), 2), test_x, test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for at least one NA in X", {
  expect_error(compute_splines_est(c(1, NA, 3), test_x, test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonnumeric x", {
  expect_error(compute_splines_est(test_X, c(1, 'a', 3), test_estCov, test_p, test_m, maxLag = 2))
  expect_error(compute_splines_est(test_X, c(1, 1i, 3), test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonvector x", {
  expect_error(compute_splines_est(test_X, matrix(c(1, 2, 3, 4), 2), test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for at least one NA in x", {
  expect_error(compute_splines_est(test_X, c(1, NA, 3), test_estCov, test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonnumeric maxLag", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, test_m, maxLag = 'a'))
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, test_m, maxLag = 1i))
})

test_that("compute_splines_est() fails for maxLag < 0", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, test_m, maxLag =-1))
})

test_that("compute_splines_est() fails for maxLag >= length(X)", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, test_m, maxLag = length(test_X)))
})

test_that("compute_splines_est() fails for nonnumeric estCov", {
  expect_error(compute_splines_est(test_X, test_x, c(1, 'a', 3), test_p, test_m, maxLag = 2))
  expect_error(compute_splines_est(test_X, test_x, c(1, 1i, 3), test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonvector estCov", {
  expect_error(compute_splines_est(test_X, test_x, matrix(c(1, 2, 3, 4), 2), test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for at least one NA in estCov", {
  expect_error(compute_splines_est(test_X, test_x, c(1, NA, 3), test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails when length(estCov) != maxLag", {
  expect_error(compute_splines_est(test_X, test_x, c(1, 2, 3, 4), test_p, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonnumeric p", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, 'a', test_m, maxLag = 2))
  expect_error(compute_splines_est(test_X, test_x, test_estCov, 1i, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for p < 0", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, -1, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for noninteger p", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, 1.1, test_m, maxLag = 2))
})

test_that("compute_splines_est() fails for nonnumeric m", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, 'a', maxLag = 2))
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, 1i, maxLag = 2))
})

test_that("compute_splines_est() fails for m < 0", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, -1, maxLag = 2))
})

test_that("compute_splines_est() fails for noninteger m", {
  expect_error(compute_splines_est(test_X, test_x, Test_estCov, test_p, 1.1, maxLag = 2))
})

test_that("compute_splines_est() fails for length(initial_pars) !=  m + p", {
  expect_error(compute_splines_est(test_X, test_x, test_estCov, test_p, test_m, initial_pars = c(), maxLag = 2))
})
