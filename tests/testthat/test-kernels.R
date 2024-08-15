# Gaussian kernel
test_that("kernel_gaussian() works", {
  expect_equal(kernel_gaussian(c(0, 1, Inf), 1), c(1, exp(-(1^2)), 0))
})

test_that("kernel_gaussian() fails for empty x, theta=1", {
  expect_error(kernel_gaussian(c(), 1))
})

test_that("kernel_gaussian() fails for theta <= 0", {
  expect_error(kernel_gaussian(c(0, 1, Inf), 0))
  expect_error(kernel_gaussian(c(0, 1, Inf), -0.1))
})

test_that("kernel_gaussian() fails for any negative x, theta=1", {
  expect_error(kernel_gaussian(c(0, 1, -1, -2), 1))
})

# Exponential kernel
test_that("kernel_exponential() works", {
  expect_equal(kernel_exponential(c(0, 1, Inf), 1), c(1, exp(-1), 0))
})

test_that("kernel_exponential() fails for empty x, theta=1", {
  expect_error(kernel_exponential(c(), 1))
})

test_that("kernel_exponential() fails for theta <= 0", {
  expect_error(kernel_exponential(c(0, 1, Inf), 0))
  expect_error(kernel_exponential(c(0, 1, Inf), -0.1))
})

test_that("kernel_exponential() fails for any negative x, theta=1", {
  expect_error(kernel_exponential(c(0, 1, -1, -2), 1))
})

# Wave kernel
test_that("kernel_wave() works", {
  expect_equal(kernel_wave(c(0, 1, Inf), 1), c(1, (1/1) * sin(1 / 1), 0))
})

test_that("kernel_wave() fails for empty x, theta=1", {
  expect_error(kernel_wave(c(), 1))
})

test_that("kernel_wave() fails for theta <= 0", {
  expect_error(kernel_wave(c(0, 1, Inf), 0))
  expect_error(kernel_wave(c(0, 1, Inf), -0.1))
})

test_that("kernel_wave() fails for any negative x, theta=1", {
  expect_error(kernel_wave(c(0, 1, -1, -2), 1))
})

# Rational Quadratic kernel
test_that("kernel_rational_quadratic() works", {
  expect_equal(kernel_rational_quadratic(c(0, 1, Inf), 1), c((1 - (0^2 / (0^2 + 1))), (1 - (1^2 / (1^2 + 1))), 0))
})

test_that("kernel_rational_quadratic() fails for empty x, theta=1", {
  expect_error(kernel_rational_quadratic(c(), 1))
})

test_that("kernel_rational_quadratic() fails for theta <= 0", {
  expect_error(kernel_rational_quadratic(c(0, 1, Inf), 0))
  expect_error(kernel_rational_quadratic(c(0, 1, Inf), -0.1))
})

test_that("kernel_rational_quadratic() fails for any negative x, theta=1", {
  expect_error(kernel_rational_quadratic(c(0, 1, -1, -2), 1))
})

# Spherical kernel
test_that("kernel_spherical() works", {
  expect_equal(kernel_spherical(c(0, 1, Inf), 2), c(1, 1 - ((3/2) * (1 / 2)) + ((1/2) * (1 / 2)^3), 0))
})

test_that("kernel_spherical() fails for empty x, theta=1", {
  expect_error(kernel_spherical(c(), 1))
})

test_that("kernel_spherical() fails for theta <= 0", {
  expect_error(kernel_spherical(c(0, 1, Inf), 0))
  expect_error(kernel_spherical(c(0, 1, Inf), -0.1))
})

test_that("kernel_spherical() fails for any negative x, theta=1", {
  expect_error(kernel_spherical(c(0, 1, -1, -2), 1))
})

# Circular kernel
test_that("kernel_circular() works", {
  expect_equal(kernel_circular(c(0, 1, Inf), 2), c(1, ((2 / pi) * acos(1 / 2)) - ((2 / pi) * (1 / 2) * sqrt(1 - (1 / 2)^2)), 0))
})

test_that("kernel_circular() fails for empty x, theta=1", {
  expect_error(kernel_circular(c(), 1))
})

test_that("kernel_circular() fails for theta <= 0", {
  expect_error(kernel_circular(c(0, 1, Inf), 0))
  expect_error(kernel_circular(c(0, 1, Inf), -0.1))
})

test_that("kernel_circular() fails for any negative x, theta=1", {
  expect_error(kernel_circular(c(0, 1, -1, -2), 1))
})

# Bessel kernel
test_that("kernel_bessel_j() works", {
  expect_equal(kernel_bessel_j(c(0, 1, Inf), 1, 0, 2), c(1, (2^0) * gamma(0 + 1) * (besselJ(1 / 1, 0) / ((1 / 1)^0)), 0))
})

test_that("kernel_bessel_j() fails for empty x, theta=1", {
  expect_error(kernel_bessel_j(c(), 1, 0, 2))
})

test_that("kernel_bessel_j() fails for theta <= 0", {
  expect_error(kernel_bessel_j(c(0, 1, Inf), 0, 0, 2))
  expect_error(kernel_bessel_j(c(0, 1, Inf), -0.1, 0, 2))
})

test_that("kernel_bessel_j() fails for any negative x, theta=1", {
  expect_error(kernel_bessel_j(c(0, 1, -1, -2), 1, 0, 2))
})

test_that("kernel_bessel_j() fails for nu <= (dim / 2) - 1, theta=1, dim=2", {
  expect_error(kernel_bessel_j(c(0, 1, Inf), 1, -0.5, 2))
})

# Matérn kernel
test_that("kernel_matern() works", {
  expect_equal(kernel_matern(c(0, 1, Inf), 1, 1), c(1, (((sqrt(2 * 1) * 1 / 1)^1) / (2^(1 - 1) * gamma(1))) * besselK(sqrt(2 * 1) * 1 / 1, 1), 0))
})

test_that("kernel_matern() fails for empty x, theta=1", {
  expect_error(kernel_matern(c(), 1, 1))
})

test_that("kernel_matern() fails for theta <= 0", {
  expect_error(kernel_matern(c(0, 1, Inf), 0, 1))
  expect_error(kernel_matern(c(0, 1, Inf), -0.1, 1))
})

test_that("kernel_matern() fails for any negative x, theta=1", {
  expect_error(kernel_matern(c(0, 1, -1, -2), 1, 0))
})

test_that("kernel_matern() fails for nu <= 0, theta=1", {
  expect_error(kernel_matern(c(0, 1, Inf), 1, 0))
  expect_error(kernel_matern(c(0, 1, Inf), 1, -0.5))
})

# Cauchy kernel
test_that("kernel_cauchy() works", {
  expect_equal(kernel_cauchy(c(0, 1, Inf), 1, 1, 1), c(1, 0.5, 0))
})

test_that("kernel_cauchy() fails for empty x, theta=1", {
  expect_error(kernel_cauchy(c(), 1, 1, 1))
})

test_that("kernel_cauchy() fails for theta < 0", {
  expect_error(kernel_cauchy(c(0, 1, Inf), 0, 1, 1))
  expect_error(kernel_cauchy(c(0, 1, Inf), -0.1, 1, 1))
})

test_that("kernel_cauchy() fails for any negative x, theta=1", {
  expect_error(kernel_cauchy(c(0, 1, -1, -2), 1, 0))
})

test_that("kernel_cauchy() fails for alpha <= 0, theta=1, beta=1", {
  expect_error(kernel_cauchy(c(0, 1, Inf), 1, 0, 1))
  expect_error(kernel_cauchy(c(0, 1, Inf), 1, -0.1, 1))
})

test_that("kernel_cauchy() fails for alpha > 2, theta=1, beta=1", {
  expect_error(kernel_cauchy(c(0, 1, Inf), 1, 3, 1))
})

test_that("kernel_cauchy() fails for beta < 0, theta=1, alpha=1", {
  expect_error(kernel_cauchy(c(0, 1, Inf), 1, 1, -0.1))
})
