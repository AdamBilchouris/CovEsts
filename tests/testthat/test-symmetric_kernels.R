# Gaussian kernel
test_that("kernel_symm_ec(, \"gaussian\", ) works", {
  expect_equal(kernel_symm_ec(c(0, 1, Inf), "gaussian", c(1)), c(sqrt(pi)^(-1)*1, sqrt(pi)^(-1)*exp(-(1^2)), 0))
})

test_that("kernel_symm_ec(, \"gaussian\", ) fails for empty x, theta=1", {
  expect_error(kernel_symm_ec(c(), "gaussian", c(1)))
})

test_that("kernel_gaussian(, \"gaussian\", ) fails for theta <= 0", {
  expect_error(kernel_symm_ec(c(0, 1, Inf), "gaussian", c(0)))
  expect_error(kernel_symm_ec(c(0, 1, Inf), "gaussian", c(-0.1)))
})

test_that("kernel_symm_ec(, \"gaussian\", ) fails for any NA x, theta=1", {
  expect_error(kernel_symm_ec(c(0, 1, NA), "gaussian", c(1)))
})

# Wave kernel
test_that("kernel_symm_ec(, \"wave\", ) works", {
  expect_equal(kernel_symm_ec(c(0, 1, Inf), "wave", c(1)), c((sqrt(1^(-2)) / pi) * 1, (sqrt(1^(-2)) / pi) * (1/1) * sin(1 / 1), 0))
})

test_that("kernel_symm_ec(, \"wave\", ) fails for empty x, theta=1", {
  expect_error(kernel_symm_ec(c(), "wave", c(1)))
})

test_that("kernel_symm_ec(, \"wave\", ) fails for theta <= 0", {
  expect_error(kernel_symm_ec(c(0, 1, Inf), "wave", c(0)))
  expect_error(kernel_symm_ec(c(0, 1, Inf), "wave", c(-0.1)))
})

test_that("kernel_symm_ec(, \"wave\", ) fails for any NA x, theta=1", {
  expect_error(kernel_symm_ec(c(0, 1, NA), "wave", c(1)))
})

# Rational Quadratic kernel
test_that("kernel_symm_ec(, \"rational_quadratic\", ) works", {
  expect_equal(kernel_symm_ec(c(0, 1, Inf), "rational_quadratic", c(1)), c(((pi * sqrt(1))^(-1)) * (1 - (0^2 / (0^2 + 1))), ((pi * sqrt(1))^(-1)) * (1 - (1^2 / (1^2 + 1))), 0))
})

test_that("kernel_symm_ec(, \"rational_quadratic\", ) fails for empty x, theta=1", {
  expect_error(kernel_symm_ec(c(), "rational_quadratic", c(1)))
})

test_that("kernel_symm_ec(, \"rational_quadratic\", ) fails for theta <= 0", {
  expect_error(kernel_symm_ec(c(0, 1, Inf), "rational_quadratic", c(0)))
  expect_error(kernel_symm_ec(c(0, 1, Inf), "rational_quadratic",  c(-0.1)))
})

test_that("kernel_symm_ec(, \"rational_quadratic\", ) fails for any NA x, theta=1", {
  expect_error(kernel_symm_ec(c(0, 1, NA), "rational_quadratic", c(1)))
})

# Bessel kernel
test_that("kernel_symm_ec(, \"bessel_j\", ) works", {
  expect_equal(kernel_symm_ec(c(0, 1, Inf), "bessel_j", c(1, 0, 2)), c((gamma((1/2) + 0) / (2*sqrt(pi) * 1 * gamma(1 + 0))) * 1,
              (gamma((1/2) + 0) / (2*sqrt(pi) * 1 * gamma(1 + 0))) * ((2^0) * gamma(0 + 1) * (besselJ(1 / 1, 0) / ((1 / 1)^0))), 0))
})

test_that("kernel_symm_ec(, \"bessel_j\", ) fails for empty x, theta=1", {
  expect_error(kernel_symm_ec(c(), "bessel_j", c(1, 0, 2)))
})

test_that("kernel_symm_ec(, \"bessel_j\", ) fails for theta <= 0", {
  expect_error(kernel_symm_ec(c(0, 1, Inf), "bessel_j", c(0, 0, 2)))
  expect_error(kernel_symm_ec(c(0, 1, Inf), "bessel_j", c(-0.1, 0, 2)))
})

test_that("kernel_symm_ec(, \"bessel_j\", ) fails for any NA x, theta=1", {
  expect_error(kernel_symm_ec(c(0, 1, NA), "bessel_j", c(1)))
})

test_that("kernel_symm_ec(, \"gaussian\") works", {
  expect_equal(kernel_symm_ec(c(0, 1, 2), "gaussian"), c(0.56418958, 0.20755375, 0.01033349))
})

test_that("kernel_symm_ec(, \"gaussian\") fails for at least one NA in x", {
  expect_error(kernel_symm_ec(c(0, NA, 2), "gaussian"))
})

test_that("kernel_symm_ec(, \"gaussian\") fails for nonnumeric x", {
  expect_error(kernel_symm_ec(c(0, 'a', 2), "gaussian"))
  expect_error(kernel_symm_ec(c(0, 1i, 2), "gaussian"))
})

test_that("kernel_symm_ec(, \"gaussian\") fails empty x", {
  expect_error(kernel_symm_ec(c(), "gaussian"))
})
