# Gaussian kernel
test_that("kernel(, \"gaussian\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "gaussian", c(1)), c(1, exp(-(1^2)), 0))
})

test_that("kernel(, \"gaussian\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "gaussian", c(1)))
})

test_that("kernel(, \"gaussian\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "gaussian", c(0)))
  expect_error(kernel(c(0, 1, Inf), "gaussian", c(-0.1)))
})

test_that("kernel(, \"exponential\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "gaussian", c(1)))
})

# Exponential kernel
test_that("kernel(, \"exponential\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "exponential", c(1)), c(1, exp(-1), 0))
})

test_that("kernel(, \"exponential\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "exponential", c(1)))
})

test_that("kernel(, \"exponential\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "exponential", c(0)))
  expect_error(kernel(c(0, 1, Inf), "exponential", c(-0.1)))
})

test_that("kernel(, \"exponential\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "exponential", c(1)))
})

# Wave kernel
test_that("kernel(, \"wave\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "wave", c(1)), c(1, (1/1) * sin(1 / 1), 0))
})

test_that("kernel(, \"wave\", )  fails for empty x, theta=1", {
  expect_error(kernel(c(), "wave", c(1)))
})

test_that("kernel(, \"wave\", )  fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "wave", c(0)))
  expect_error(kernel(c(0, 1, Inf), "wave", c(-0.1)))
})

test_that("kernel(, \"wave\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "wave", c(1)))
})

# Rational Quadratic kernel
test_that("kernel(, \"rational_quadratic\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "rational_quadratic", c(1)), c((1 - (0^2 / (0^2 + 1))), (1 - (1^2 / (1^2 + 1))), 0))
})

test_that("kernel(, \"rational_quadratic\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "rational_quadratic", c(1)))
})

test_that("kernel(, \"rational_quadratic\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "rational_quadratic", c(0)))
  expect_error(kernel(c(0, 1, Inf), "rational_quadratic", (-0.1)))
})

test_that("kernel(, \"rational_quadratic\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "rational_quadratic", c(1)))
})

# Spherical kernel
test_that("kernel(, \"spherical\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "spherical", c(2)), c(1, 1 - ((3/2) * (1 / 2)) + ((1/2) * (1 / 2)^3), 0))
})

test_that("kernel(, \"spherical\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "spherical", c(1)))
})

test_that("kernel(, \"spherical\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "spherical", c(0)))
  expect_error(kernel(c(0, 1, Inf), "spherical", c(-0.1)))
})

test_that("kernel(, \"spherical\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "spherical", c(1)))
})

# Circular kernel
test_that("kernel(, \"circular\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "circular", c(2)), c(1, ((2 / pi) * acos(1 / 2)) - ((2 / pi) * (1 / 2) * sqrt(1 - (1 / 2)^2)), 0))
})

test_that("kernel(, \"circular\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "circular", c(1)))
})

test_that("kernel(, \"circular\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "circular", c(0)))
  expect_error(kernel(c(0, 1, Inf), "circular", c(-0.1)))
})

test_that("kernel(, \"circular\", )fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "circular", c(1)))
})

# Bessel kernel
test_that("kernel(, \"bessel_j\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "bessel_j", c(1, 0, 2)), c(1, (2^0) * gamma(0 + 1) * (besselJ(1 / 1, 0) / ((1 / 1)^0)), 0))
})

test_that("kernel(, \"bessel_j\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), 1, 0, 2))
})

test_that("kernel(, \"bessel_j\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "bessel_j",  c(0, 0, 2)))
  expect_error(kernel(c(0, 1, Inf), "bessel_j",  c(-0.1, 0, 2)))
})

test_that("kernel(, \"bessel_j\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "bessel_j",  c(1, 0, 2)))
})

test_that("kernel(, \"bessel_j\", ) fails for nu <= (dim / 2) - 1, theta=1, dim=2", {
  expect_error(kernel(c(0, 1, Inf), "bessel_j", c(1, -0.5, 2)))
})

# Matérn kernel
test_that("kernel(, \"matern\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "matern", c(1, 1)), c(1, (((sqrt(2 * 1) * 1 / 1)^1) / (2^(1 - 1) * gamma(1))) * besselK(sqrt(2 * 1) * 1 / 1, 1), 0))
})

test_that("kernel(, \"matern\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "matern", c(1, 1)))
})

test_that("kernel(, \"matern\", ) fails for theta <= 0", {
  expect_error(kernel(c(0, 1, Inf), "matern", c(0, 1)))
  expect_error(kernel(c(0, 1, Inf), "matern", c(-0.1, 1)))
})

test_that("kernel(, \"matern\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "matern", c(1, 0)))
})

test_that("kernel(, \"matern\", ) fails for nu <= 0, theta=1", {
  expect_error(kernel(c(0, 1, Inf), "matern", c(1, 0)))
  expect_error(kernel(c(0, 1, Inf), "matern", c(1, -0.5)))
})

# Cauchy kernel
test_that("kernel(, \"cauchy\", ) works", {
  expect_equal(kernel(c(0, 1, Inf), "cauchy", c(1, 1, 1)), c(1, 0.5, 0))
})

test_that("kernel(, \"cauchy\", ) fails for empty x, theta=1", {
  expect_error(kernel(c(), "cauchy", c(1, 1, 1)))
})

test_that("kernel(, \"cauchy\", ) fails for theta < 0", {
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(0, 1, 1)))
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(-0.1, 1, 1)))
})

test_that("kernel(, \"cauchy\", ) fails for any negative x, theta=1", {
  expect_error(kernel(c(0, 1, -1, -2), "cauchy", c(1, 1, 1)))
})

test_that("kernel(, \"cauchy\", ) fails for alpha <= 0, theta=1, beta=1", {
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(1, 0, 1)))
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(1, -0.1, 1)))
})

test_that("kernel(, \"cauchy\", ) fails for alpha > 2, theta=1, beta=1", {
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(1, 3, 1)))
})

test_that("kernel(, \"cauchy\", ) fails for beta < 0, theta=1, alpha=1", {
  expect_error(kernel(c(0, 1, Inf), "cauchy", c(1, 1, -0.1)))
})
