# Tukey window
test_that("window_ec(, \"tukey\") works", {
  expect_equal(window_ec(c(1, 2, 3) / 3, "tukey"), c(0.25, 0.75, 1))
})

test_that("window_ec(, \"tukey\") fails for empty x", {
  expect_error(window_ec(c(), "tukey"))
})

test_that("window_ec(, \"tukey\") fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "tukey"))
})

test_that("window_ec(, \"tukey\") fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "tukey"))
})

# Triangular window
test_that("window_ec(, \"trianguar\") works", {
  expect_equal(window_ec(c(1, 2, 3) / 4, "triangular"), c(0.25, 0.5, 0.75))
})

test_that("window_ec(, \"trianguar\") fails for empty x", {
  expect_error(window_ec(c(), "triangular"))
})

test_that("window_ec(, \"trianguar\") fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "triangular"))
})

test_that("window_ec(, \"trianguar\") fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "triangular"))
})

# Sine window
test_that("window_ec(, \"sine\") works", {
  expect_equal(window_ec(c(0, 1/3, 1/2, 2/3), "sine"), c(0, sin(pi/6), sin(pi/4), sin(pi/3)))
})

test_that("window_ec(, \"sine\") fails for empty x", {
  expect_error(window_ec(c(), "sine"))
})

test_that("window_ec(, \"sine\") fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "sine"))
})

test_that("window_ec(, \"sine\") fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "sine"))
})

# Power sine window
test_that("window_ec(, \"power_sine\", ) works", {
  expect_equal(window_ec(c(0, 1/3, 1/2, 2/3), "power_sine", c(2)), c(0, 0.25, 0.5, 0.75))
})

test_that("window_ec(, \"power_sine\", ) fails for empty x", {
  expect_error(window_ec(c(), "power_sine", c(2)))
})

test_that("window_ec(, \"power_sine\", ) fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "power_sine", c(2)))
})

test_that("window_ec(, \"power_sine\", ) fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "power_sine", c(2)))
})

test_that("window_ec(, \"power_sine\", ) fails for a <= 0", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "power_sine", c(0)))
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "power_sine", c(-1)))
})

# Blackman window
test_that("window_ec(, \"blackman\", ) works", {
  expect_equal(window_ec(c(0, 1/3, 1/2, 2/3), "blackman", c(0.16)), c(-1.387779e-17, 1.300000e-01, 3.400000e-01, 6.300000e-01))
})

test_that("window_ec(, \"blackman\", ) fails for empty x", {
  expect_error(window_ec(c(), "blackman", c(2)))
})

test_that("window_ec(, \"blackman\", ) fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "blackman", c(2)))
})

test_that("window_ec(, \"blackman\", ) fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "blackman", c(2)))
})

test_that("window_ec(, \"blackman\", ) fails for a nonreal alpha.", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "blackman", c(1i)))
})

# Hann-Poisson window
test_that("window_ec(, \"hann_poisson\", ) works", {
  expect_equal(window_ec(c(0, 1/3, 1/2, 2/3), "hann_poisson", c(0.5)), c(0, (1 / 2) * (1 - cos(pi * (1/3))) * exp(- (0.5 * abs(1 - (1/3)))),
                                                                (1 / 2) * (1 - cos(pi * (1/2))) * exp(- (0.5 * abs(1 - (1/2)))),
                                                                ((1 / 2) * (1 - cos(pi * (2/3))) * exp(- (0.5 * abs(1 - (2/3)))))))
})

test_that("window_ec(, \"hann_poisson\", ) fails for empty x", {
  expect_error(window_ec(c(), "hann_poisson", c(2)))
})

test_that("window_ec(, \"hann_poisson\", ) fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "hann_poisson", c(2)))
})

test_that("window_ec(, \"hann_poisson\", ) fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "hann_poisson", c(2)))
})

test_that("window_ec(, \"blackman\", ) fails for a nonreal alpha.", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "hann_poisson", c(1i)))
})

# Welch window
test_that("window_ec(, \"welch\") works", {
  expect_equal(window_ec(c(0, 1/3, 1/2, 2/3), "welch"), c(0, 1 - ((1/3)-1)^2, 0.75, 1 - ((2/3)-1)^2))
})

test_that("window_ec(, \"welch\") fails for empty x", {
  expect_error(window_ec(c(), "welch"))
})

test_that("window_ec(, \"welch\") fails for a single negative x", {
  expect_error(window_ec(c(0.1 , 0.3, -0.02), "welch"))
})

test_that("window_ec(, \"welch\") fails for a single x above 1", {
  expect_error(window_ec(c(0.1 , 0.3, 1.02), "welch"))
})

test_that("window_symm_ec(, \"tukey\") works", {
  expect_equal(window_symm_ec(c(1, 2, 3) / 3, "tukey"), c(0.75, 0.25, 0))
})

test_that("window_symm_ec(, \"tukey\") fails for at least one NA in x", {
  expect_error(window_symm_ec(c(1, NA, 3) / 3, "tukey"))
})

test_that("window_symm_ec(, \"tukey\") fails for nonnumeric x", {
  expect_error(window_symm_ec(c(1, 'a', 3), "tukey"))
  expect_error(window_symm_ec(c(1, 1i, 3) / 3, "tukey"))
})

test_that("window_symm_ec(, \"tukey\") fails empty x", {
  expect_error(window_symm_ec(c(), "tukey"))
})

test_that('.validate_window_params() works', {
  expect_equal(.validate_window_params("tukey", c(1)), TRUE)
})

test_that('.validate_window_params() fails for name not in "tukey", "triangular", "sine", "power_sine", "blackman", "hann_poisson", "welch"', {
  expect_error(.validate_window_params("test", c(1)))
})
