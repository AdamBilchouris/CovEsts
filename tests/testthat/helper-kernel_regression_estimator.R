# This function sets variables such use in parameters.
test_x <- c(1, 2, 3)
test_X <- c(1, 2, 3)
test_meanX <- mean(test_X)
test_xij <- Xij_mat(test_X)
test_rhoT1 <- rho_T1(test_x, test_meanX, 2, 0.01, test_xij)
