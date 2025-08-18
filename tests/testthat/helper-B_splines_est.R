test_p <- 2
test_m <- 3
test_kVec <- generate_knots(test_m)
test_taus <- get_taus(test_p, test_m)
test_x <- c(1, 2, 3)
test_X <- c(1, 2, 3)

test_estCov <- standard_est(test_X, maxLag=length(test_X) - 1)
