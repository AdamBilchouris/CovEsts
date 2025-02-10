test_X <- c(1, 2, 3)
test_meanX <- mean(test_X)
test_h2n <- H2n(length(test_X), 1, "tukey")
test_n <- 3
test_taperVals_t <- taper(((1:test_n) - 1/2)/test_n, 1, "tukey")
test_taperVals_h <- taper((((1:(test_n-1)) - 1/2) +1) / test_n, 1, "tukey")
