# starting_locs
test_that("starting_locs fails for nonnumeric N", {
  expect_error(starting_locs('a', 2, 2))
  expect_error(starting_locs(1i, 2, 2))
})

test_that("starting_locs fails for N < 1", {
  expect_error(starting_locs(0, 2, 2))
})

test_that("starting_locs fails for NA N", {
  expect_error(starting_locs(NA, 2, 2))
})

test_that("starting_locs fails for nonnumeric l", {
  expect_error(starting_locs(4, 'a', 2))
  expect_error(starting_locs(4, 1i, 2))
})

test_that("starting_locs fails for l < 1", {
  expect_error(starting_locs(4, 0, 2))
})

test_that("starting_locs fails for NA l", {
  expect_error(starting_locs(4, NA, 2))
})

test_that("starting_locs fails for nonnumeric k", {
  expect_error(starting_locs(4, 2, 'a'))
  expect_error(starting_locs(4, 2, 1i))
})

test_that("starting_locs fails for k < 1", {
  expect_error(starting_locs(4, 2, 0))
})

test_that("starting_locs fails for NA k", {
  expect_error(starting_locs(4, 2, NA))
})

test_that("starting_locs fails for boot_type being neither 'moving' or 'circular'", {
  expect_error(starting_locs(4, 2, 2, 'stationary'))
})

#  bootstrap_samples
test_that("bootstrap_samples fails for nonnumeric X", {
  expect_error(bootstrap_samples(c(1,'a', 3, 4), 2, 2))
  expect_error(bootstrap_samples(c(1, 2i, 3, 4), 2, 2))
})

test_that("bootstrap_samples fails for X less than length 1", {
  expect_error(bootstrap_samples(c(), 2, 2))
})

test_that("bootstrap_samples fails for at least one NA in X", {
  expect_error(bootstrap_samples(c(1, NA, 3, 4), 2, 2))
})

test_that("bootstrap_samples fails for nonnumeric l", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 'a', 2))
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 1i, 2))
})

test_that("bootstrap_samples fails for l < 1", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 0, 2))
})

test_that("bootstrap_samples fails for NA l", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), NA, 2))
})

test_that("bootstrap_samples fails for nonnumeric k", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 2, 'a'))
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 2, 1i))
})

test_that("bootstrap_samples fails for k < 1", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 2, 0))
})

test_that("bootstrap_samples fails for NA k", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 2, NA))
})

test_that("bootstrap_samples fails for boot_type being neither 'moving' or 'circular'", {
  expect_error(bootstrap_samples(c(1, 2, 3, 4), 2, 2, 'stationary'))
})

# block_bootstrap
test_that("block_bootstrap fails for nonnumeric X", {
  expect_error(block_bootstrap(c(1,'a', 3, 4), 2))
  expect_error(block_bootstrap(c(1, 2i, 3, 4), 2))
})

test_that("block_bootstrap fails for X less than length 1", {
  expect_error(block_bootstrap(c(), 2))
})

test_that("block_bootstrap fails for at least one NA in X", {
  expect_error(block_bootstrap(c(1, NA, 3, 4), 2))
})

test_that("block_bootstrap fails for nonnumeric maxLag", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 'a'))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 1i))
})

test_that("block_bootstrap fails for maxLag not being length 1", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), c(1, 2)))
})

test_that("block_bootstrap fails for maxLag not being greater than 0", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 0))
})

test_that("block_bootstrap fails for maxLag not being less than length(X)-1", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 4))
})

test_that("block_bootstrap fails for maxLag not being an integer", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2.5))
})

test_that("block_bootstrap fails for nonnumeric x", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, x=c(1, 'a', 3, 4)))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, x=c(1, 2i, 3, 4)))
})

test_that("block_bootstrap fails for length(x) != length(X)", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, x=1))
})

test_that("block_bootstrap fails for nonnumeric n_bootstrap", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, n_bootstrap='a'))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, n_bootstrap=1i))
})

test_that("block_bootstrap fails n_bootstrap being less than or equal to zero", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, n_bootstrap=0))
})

test_that("block_bootstrap fails noninteger n_bootstrap", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, n_bootstrap=50.5))
})

test_that("block_bootstrap fails for nonnumeric l", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l='a'))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l=1i))
})

test_that("block_bootstrap fails for lenght of l not being 1", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l=c(1, 2)))
})

test_that("block_bootstrap fails for l being less than or equal to 0", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l=0))
})

test_that("block_bootstrap fails for l being greater than length(X)", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l=5))
})

test_that("block_bootstrap fails noninteger l", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, l=2.5))
})

test_that("block_bootstrap fails for nonnumeric alpha", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, alpha='a'))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, alpha=1i))
})

test_that("block_bootstrap fails for alpha being greater than 1", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, alpha=1.01))
})

test_that("block_bootstrap fails for alpha being less than 0", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, alpha=-0.01))
})

test_that("block_bootstrap fails for nonboolean plot", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, plot = 1))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, plot = 'TRUE'))
})

test_that("block_bootstrap fails for nonboolean boot_mat", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, boot_mat = 1))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, boot_mat = 'TRUE'))
})

test_that("block_bootstrap fails for nonnumeric ylim", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, ylim = c('a', 1)))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, ylim = c(1i, 1)))
})

test_that("block_bootstrap fails for ylim not of length 2", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, ylim = c()))
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, ylim = c(-1, 0, 1)))
})

test_that("block_bootstrap fails for at least type being neither autocovariance or autocorrelation", {
  expect_error(block_bootstrap(c(1, 2, 3, 4), 2, type = "covariance"))
})
