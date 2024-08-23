# Gaussian kernel
test_that("kernel_symm_gaussian() works", {
  expect_equal(kernel_symm_gaussian(c(0, 1, Inf), 1), c(sqrt(pi)^(-1)*1, sqrt(pi)^(-1)*exp(-(1^2)), 0))
})

test_that("kernel_symm_gaussian() fails for empty x, theta=1", {
  expect_error(kernel_symm_gaussian(c(), 1))
})

test_that("kernel_gaussian() fails for theta <= 0", {
  expect_error(kernel_symm_gaussian(c(0, 1, Inf), 0))
  expect_error(kernel_symm_gaussian(c(0, 1, Inf), -0.1))
})

test_that("kernel_symm_gaussian() fails for any NA x, theta=1", {
  expect_error(kernel_symm_gaussian(c(0, 1, NA), 1))
})

# Wave kernel
test_that("kernel_symm_wave() works", {
  expect_equal(kernel_symm_wave(c(0, 1, Inf), 1), c((sqrt(1^(-2)) / pi) * 1, (sqrt(1^(-2)) / pi) * (1/1) * sin(1 / 1), 0))
})

test_that("kernel_symm_wave() fails for empty x, theta=1", {
  expect_error(kernel_symm_wave(c(), 1))
})

test_that("kernel_symm_wave() fails for theta <= 0", {
  expect_error(kernel_symm_wave(c(0, 1, Inf), 0))
  expect_error(kernel_symm_wave(c(0, 1, Inf), -0.1))
})

test_that("kernel_symm_wave() fails for any NA x, theta=1", {
  expect_error(kernel_symm_wave(c(0, 1, NA), 1))
})

# Rational Quadratic kernel
test_that("kernel_symm_rational_quadratic() works", {
  expect_equal(kernel_symm_rational_quadratic(c(0, 1, Inf), 1), c(((pi * sqrt(1))^(-1)) * (1 - (0^2 / (0^2 + 1))), ((pi * sqrt(1))^(-1)) * (1 - (1^2 / (1^2 + 1))), 0))
})

test_that("kernel_symm_rational_quadratic() fails for empty x, theta=1", {
  expect_error(kernel_symm_rational_quadratic(c(), 1))
})

test_that("kernel_symm_rational_quadratic() fails for theta <= 0", {
  expect_error(kernel_symm_rational_quadratic(c(0, 1, Inf), 0))
  expect_error(kernel_symm_rational_quadratic(c(0, 1, Inf), -0.1))
})

test_that("kernel_symm_rational_quadratic() fails for any NA x, theta=1", {
  expect_error(kernel_symm_rational_quadratic(c(0, 1, NA), 1))
})


# Bessel kernel
test_that("kernel_symm_bessel_j() works", {
  expect_equal(kernel_symm_bessel_j(c(0, 1, Inf), 1, 0, 2), c((gamma((1/2) + 0) / (2*sqrt(pi) * 1 * gamma(1 + 0))) * 1,
              (gamma((1/2) + 0) / (2*sqrt(pi) * 1 * gamma(1 + 0))) * ((2^0) * gamma(0 + 1) * (besselJ(1 / 1, 0) / ((1 / 1)^0))), 0))
})

test_that("kernel_symm_bessel_j() fails for empty x, theta=1", {
  expect_error(kernel_symm_bessel_j(c(), 1, 0, 2))
})

test_that("kernel_symm_bessel_j() fails for theta <= 0", {
  expect_error(kernel_symm_bessel_j(c(0, 1, Inf), 0, 0, 2))
  expect_error(kernel_symm_bessel_j(c(0, 1, Inf), -0.1, 0, 2))
})

test_that("kernel_symm_bessel_j() fails for any NA x, theta=1", {
  expect_error(kernel_symm_bessel_j(c(0, 1, NA), 1))
})
