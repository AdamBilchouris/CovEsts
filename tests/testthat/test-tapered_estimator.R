# H2n
test_that("H2n(, , \"tukey\") works", {
  expect_equal(H2n(3, 1, "tukey"), 1.5)
})

test_that("H2n(, , \"tukey\") fails for nonnumeric n", {
  expect_error(H2n('a', 1, "tukey"))
  expect_error(H2n(1i, 1, "tukey"))
})

test_that("H2n(, , \"tukey\") fails for n < 1", {
  expect_error(H2n(0, 1, "tukey"))
})

test_that("H2n(, , \"tukey\") fails for noninteger n", {
  expect_error(H2n(2.95, 1, "tukey"))
})

test_that("H2n() fails for nonnumeric rho", {
  expect_error(H2n(3, 'a', "tukey"))
  expect_error(H2n(3, 1i, "tukey"))
})

test_that("H2n() fails for negative rho", {
  expect_error(H2n(3, -0.01, "tukey"))
})

test_that("H2n() fails for rho > 1", {
  expect_error(H2n(3, 1.01, "tukey"))
})

test_that("H2n() fails for rho that is not a single number", {
  expect_error(H2n(0.5, c(1, 2), "tukey"))
})

test_that("H2n() fails for window_name that is not one of 'tukey', 'triangular', 'power_sine', 'blackman_window', 'hann_poisson', 'welch'", {
  expect_error(H2n(3, 1, "window"))
})

test_that("H2n() fails for nonboolean custom_window", {
  expect_error(H2n(3, 1, "window", custom_window = 1))
  expect_error(H2n(3, 1, "window", custom_window = 'TRUE'))
})

# taper
test_that("taper(, , \"tukey\") works", {
  expect_equal(taper(c(0.25, 0.5, 0.75), 1, "tukey"), c(0.5, 1, 0.5))
})

test_that("taper(, , \"tukey\") fails for nonnumeric x", {
  expect_error(taper(c(0.25, 'a', 0.75), 1, "tukey"))
  expect_error(taper(c(0.25, 1i, 0.75), 1, "tukey"))
})

test_that("taper(, , \"tukey\") fails for empty x", {
  expect_error(taper(c(), 1, "tukey"))
})

test_that("taper(, , \"tukey\") fails for any NAs in x", {
  expect_error(taper(c(0.25, NA, 0.75), 1, "tukey"))
})

test_that("taper() fails for nonnumeric rho", {
  expect_error(taper(3, 'a', "tukey"))
  expect_error(taper(3, 1i, "tukey"))
})

test_that("taper() fails for negative rho", {
  expect_error(taper(3, -0.01, "tukey"))
})

test_that("taper() fails for rho > 1", {
  expect_error(taper(3, 1.01, "tukey"))
})

test_that("taper() fails for rho that is not a single number", {
  expect_error(taper(0.5, c(1, 2), "tukey"))
})

test_that("taper() fails for nonboolean custom_window", {
  expect_error(taper(3, 1, "window", custom_window = 1))
  expect_error(taper(3, 1, "window", custom_window = 'TRUE'))
})

# tapered_est
test_that("tapered_est(, , , \"tukey\") works", {
  expect_equal(tapered_est(test_X, 1, "tukey", maxLag = 2)$acf, c(1/12, 0, -1/24))
})

test_that("tapered_est() fails for nonnumeric X", {
  expect_error(tapered_est(c(1, 'a', 3), 1, "tukey", maxLag = 2))
  expect_error(tapered_est(c(1, 1i, 3), 1, "tukey", maxLag = 2))
})

test_that("tapered_est() fails for empty X", {
  expect_error(tapered_est(c(), 1, "tukey", maxLag = 2))
})

test_that("tapered_est() fails for at least one NA in X", {
  expect_error(tapered_est(c(1, NA, 3), 1, "tukey", maxLag = 2))
})

test_that("tapered_est() fails for at nonnumeric maxLag", {
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = 'a'))
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = 1i))
})

test_that("tapered_est() fails for maxLag whose length is not 1", {
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = c()))
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = c(1, 2)))
})

test_that("tapered_est() fails for maxLag less than 1", {
  expect_error(tapered_est(test_X, 1, "tukey", , maxLag = 0))
})

test_that("tapered_est() fails for maxLag greater than or equal to length(X)", {
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = 3))
})

test_that("tapered_est() fails for noninteger maxLag", {
  expect_error(tapered_est(test_X, 1, "tukey", maxLag = 1.2))
})

test_that("tapered_est() fails for nonnumeric rho", {
  expect_error(tapered_est(test_X, 'a', "tukey", maxLag = 2))
  expect_error(tapered_est(test_X, 1i, "tukey", maxLag = 2))
})

test_that("tapered_est() fails for rho not of length 1", {
  expect_error(tapered_est(test_X, c(), "tukey", maxLag = 2))
  expect_error(tapered_est(test_X, c(1, 2), "tukey", maxLag = 2))
})

test_that("tapered_est() fails for rho not in (0, 1]", {
  expect_error(tapered_est(test_X, -0.01, "tukey", maxLag = 2))
  expect_error(tapered_est(test_X, 1.01, "tukey", maxLag = 2))
})

test_that("tapered_est() fails for nonboolean custom_window", {
  expect_error(tapered_est(test_X, 1, "tukey", custom_window = 1, , maxLag = 2))
  expect_error(tapered_est(test_X, 1, "tukey", custom_window = 'TRUE', maxLag = 2))
})

test_that("tapered_est() fails for nonnumeric meanX", {
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, meanX = 'a'))
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, meanX = 1i))
})

test_that("tapered_est() fails for NA meanX", {
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, meanX = NA))
})

test_that("tapered_est() fails for nonnumeric x", {
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, x = c(1, 'a', 3)))
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, x = c(1, 1i, 3)))
})

test_that("tapered_est() fails for empty x", {
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, x = c()))
})

test_that("tapered_est() fails for at least one NA in x", {
  expect_error(tapered_est(c(1, 2, 3), 1, "tukey", maxLag = 2, x = c(1, NA, 3)))
})
