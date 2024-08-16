# Tukey window
test_that("tukey_window() works", {
  expect_equal(tukey_window(c(1, 2, 3) / 3), c(0.25, 0.75, 1))
})

test_that("tukey_window() fails for empty x", {
  expect_error(tukey_window(c()))
})

test_that("tukey_window() fails for a single negative x", {
  expect_error(tukey_window(c(0.1 , 0.3, -0.02)))
})

test_that("tukey_window() fails for a single x above 1", {
  expect_error(tukey_window(c(0.1 , 0.3, 1.02)))
})

# Triangular window
test_that("triangular_window() works", {
  expect_equal(triangular_window(c(1, 2, 3) / 4), c(0.25, 0.5, 0.75))
})

test_that("triangular_window() fails for empty x", {
  expect_error(triangular_window(c()))
})

test_that("triangular_window() fails for a single negative x", {
  expect_error(triangular_window(c(0.1 , 0.3, -0.02)))
})

test_that("triangular_window() fails for a single x above 1", {
  expect_error(triangular_window(c(0.1 , 0.3, 1.02)))
})

# Sine window
test_that("sine_window() works", {
  expect_equal(sine_window(c(0, 1/3, 1/2, 2/3)), c(0, sin(pi/6), sin(pi/4), sin(pi/3)))
})

test_that("sine_window() fails for empty x", {
  expect_error(sine_window(c()))
})

test_that("sine_window() fails for a single negative x", {
  expect_error(sine_window(c(0.1 , 0.3, -0.02)))
})

test_that("sine_window() fails for a single x above 1", {
  expect_error(sine_window(c(0.1 , 0.3, 1.02)))
})

# Power sine window
test_that("power_sine_window() works", {
  expect_equal(power_sine_window(c(0, 1/3, 1/2, 2/3), 2), c(0, 0.25, 0.5, 0.75))
})

test_that("power_sine_window() fails for empty x", {
  expect_error(power_sine_window(c(), 2))
})

test_that("power_sine_window() fails for a single negative x", {
  expect_error(power_sine_window(c(0.1 , 0.3, -0.02), 2))
})

test_that("power_sine_window() fails for a single x above 1", {
  expect_error(power_sine_window(c(0.1 , 0.3, 1.02), 2))
})

# Blackman window
test_that("blackman_window() works", {
  expect_equal(blackman_window(c(0, 1/3, 1/2, 2/3), 2), c(0, -1.25, -1.5, -0.75))
})

test_that("blackman_window() fails for empty x", {
  expect_error(blackman_window(c(), 2))
})

test_that("blackman_window() fails for a single negative x", {
  expect_error(blackman_window(c(0.1 , 0.3, -0.02), 2))
})

test_that("blackman_window() fails for a single x above 1", {
  expect_error(blackman_window(c(0.1 , 0.3, 1.02), 2))
})

# Hann-Poisson window
test_that("hann_poisson_window() works", {
  expect_equal(hann_poisson_window(c(0, 1/3, 1/2, 2/3), 0.5), c(0, (1 / 2) * (1 - cos(pi * (1/3))) * exp(- (0.5 * abs(1 - (1/3)))),
                                                                (1 / 2) * (1 - cos(pi * (1/2))) * exp(- (0.5 * abs(1 - (1/2)))),
                                                                ((1 / 2) * (1 - cos(pi * (2/3))) * exp(- (0.5 * abs(1 - (2/3)))))))
})

test_that("hann_poisson_window() fails for empty x", {
  expect_error(hann_poisson_window(c(), 2))
})

test_that("hann_poisson_window() fails for a single negative x", {
  expect_error(hann_poisson_window(c(0.1 , 0.3, -0.02), 2))
})

test_that("hann_poisson_window() fails for a single x above 1", {
  expect_error(hann_poisson_window(c(0.1 , 0.3, 1.02), 2))
})

# Welch window
test_that("welch_window() works", {
  expect_equal(welch_window(c(0, 1/3, 1/2, 2/3)), c(0, 1 - ((1/3)-1)^2, 0.75, 1 - ((2/3)-1)^2))
})

test_that("welch_window() fails for empty x", {
  expect_error(welch_window(c(), 2))
})

test_that("welch_window() fails for a single negative x", {
  expect_error(welch_window(c(0.1 , 0.3, -0.02), 2))
})

test_that("welch_window() fails for a single x above 1", {
  expect_error(welch_window(c(0.1 , 0.3, 1.02), 2))
})
