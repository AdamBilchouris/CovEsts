test_p <- 2
test_m <- 3
test_kVec <- generate_knots(test_m)
test_taus <- get_all_tau(test_p, test_m)
test_x <- c(1, 2, 3)
test_X <- c(1, 2, 3)

test_estCov <- standard_est(test_X, 1)

# wrapper function that sets the seed just within this function.
# https://stackoverflow.com/a/56192474
# splines_with_seed <- function(seed, X, x, maxLag, estCov, p, m, inital_pars = c(), control=list('maxit' = 1000)) {
#   set.seed(seed)
#   return(splines_est(X, x, maxLag, estCov, p, m, inital_pars, control))
# }
