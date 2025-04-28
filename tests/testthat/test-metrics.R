# area_between
test_that("area_between() works", {
  expect_equal(area_between(test_estCov1, test_estCov2), 2.10491826, tolerance = sqrt(.Machine$double.eps))
})

test_that("area_between() works with a lag vector", {
  expect_equal(area_between(test_estCov1, test_estCov2, test_x_metrics), 0.210491826, tolerance = sqrt(.Machine$double.eps))
})

test_that("area_between() fails for nonnumeric est1", {
  expect_error(area_between(c('a', test_estCov1[-1]), test_estCov2))
  expect_error(area_between(c(1i, test_estCov1[-1]), test_estCov2))
})

test_that("area_between() fails for nonnumeric est2", {
  expect_error(area_between(test_estCov1, c('a', test_estCov2[-1])))
  expect_error(area_between(test_estCov1, c(1i, test_estCov2[-1])))
})

test_that("area_between() fails for empty est1", {
  expect_error(area_between(c(), test_estCov2))
})

test_that("area_between() fails for empty est2", {
  expect_error(area_between(test_estCov1, c()))
})

test_that("area_between() fails for at least on NA in est1", {
  expect_error(area_between(c(NA, test_estCov1[-1]), test_estCov2))
})

test_that("area_between() fails for at least on NA in est2", {
  expect_error(area_between(test_estCov1, c(NA, test_estCov2[-1])))
})

test_that("area_between() fails if est1 is of a different length", {
  expect_error(area_between(test_estCov1[-1], test_estCov2))
  expect_error(area_between(test_estCov1[-1], test_estCov2, test_x_metrics))
})

test_that("area_between() fails if est2 is of a different length", {
  expect_error(area_between(test_estCov1, test_estCov2[-1]))
  expect_error(area_between(test_estCov1, test_estCov2[-1], test_x_metrics))
})

test_that("area_between() fails if plot is not a boolean", {
  expect_error(area_between(test_estCov1, test_estCov2, plot = 1))
  expect_error(area_between(test_estCov1, test_estCov2, plot = 'TRUE'))
})

test_that("area_between() fails if lags is not a vector", {
  expect_error(area_between(test_estCov1, test_estCov2, matrix(c(1,2,3,4), 2)))
})

test_that("area_between() fails if lags is of a different length", {
  expect_error(area_between(test_estCov1, test_estCov2, test_x_metrics[-1]))
})

# max_distance
test_that("max_distance() works", {
  expect_equal(max_distance(test_estCov1, test_estCov2), 0.1818237429, tolerance = sqrt(.Machine$double.eps))
})

test_that("max_distance() fails for nonnumeric est1", {
  expect_error(max_distance(c('a', test_estCov1[-1]), test_estCov2))
  expect_error(max_distance(c(1i, test_estCov1[-1]), test_estCov2))
})

test_that("max_distance() fails for nonnumeric est2", {
  expect_error(max_distance(test_estCov1, c('a', test_estCov2[-1])))
  expect_error(max_distance(test_estCov1, c(1i, test_estCov2[-1])))
})

test_that("max_distance() fails for empty est1", {
  expect_error(max_distance(c(), test_estCov2))
})

test_that("max_distance() fails for empty est2", {
  expect_error(max_distance(test_estCov1, c()))
})

test_that("max_distance() fails for at least on NA in est1", {
  expect_error(max_distance(c(NA, test_estCov1[-1]), test_estCov2))
})

test_that("max_distance() fails for at least on NA in est2", {
  expect_error(max_distance(test_estCov1, c(NA, test_estCov2[-1])))
})

test_that("max_distance() fails if est1 is of a different length", {
  expect_error(max_distance(test_estCov1[-1], test_estCov2))
})

test_that("max_distance() fails if est2 is of a different length", {
  expect_error(max_distance(test_estCov1, test_estCov2[-1]))
})

test_that("max_distance() fails if plot is not a boolean", {
  expect_error(max_distance(test_estCov1, test_estCov2, plot = 1))
  expect_error(max_distance(test_estCov1, test_estCov2, plot = 'TRUE'))
})

# create_cyclic_matrix
test_that("create_cyclic_matrix() works", {
  expect_equal(create_cyclic_matrix(c(1, 2)), matrix(c(1, 2, 2, 1), 2))
})

test_that("create_cyclic_matrix() fails for nonnumeric v", {
  expect_error(create_cyclic_matrix(c(1, 'a', 3)))
  expect_error(create_cyclic_matrix(c(1, 1i, 3)))
})

test_that("create_cyclic_matrix() fails for at least one NA in v", {
  expect_error(create_cyclic_matrix(c(1, NA, 3)))
})


test_that("create_cyclic_matrix() fails if v is empty", {
  expect_error(create_cyclic_matrix(c()))
})

# spectral_norm
test_that("spectral_norm() works", {
  expect_equal(spectral_norm(test_estCov1, test_estCov2), 1.661283524, tolerance = sqrt(.Machine$double.eps))
})

test_that("spectral_norm() fails for nonnumeric est1", {
  expect_error(spectral_norm(c('a', test_estCov1[-1]), test_estCov2))
  expect_error(spectral_norm(c(1i, test_estCov1[-1]), test_estCov2))
})

test_that("spectral_norm() fails for nonnumeric est2", {
  expect_error(spectral_norm(test_estCov1, c('a', test_estCov2[-1])))
  expect_error(spectral_norm(test_estCov1, c(1i, test_estCov2[-1])))
})

test_that("spectral_norm() fails for empty est1", {
  expect_error(spectral_norm(c(), test_estCov2))
})

test_that("spectral_norm() fails for empty est2", {
  expect_error(spectral_norm(test_estCov1, c()))
})

test_that("spectral_norm() fails for at least on NA in est1", {
  expect_error(spectral_norm(c(NA, test_estCov1[-1]), test_estCov2))
})

test_that("spectral_norm() fails for at least on NA in est2", {
  expect_error(spectral_norm(test_estCov1, c(NA, test_estCov2[-1])))
})

test_that("spectral_norm() fails if est1 is of a different length", {
  expect_error(spectral_norm(test_estCov1[-1], test_estCov2))
})

test_that("spectral_norm() fails if est2 is of a different length", {
  expect_error(spectral_norm(test_estCov1, test_estCov2[-1]))
})
