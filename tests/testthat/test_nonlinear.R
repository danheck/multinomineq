library(testthat)


test_that("nonlinear Gibbs/BF gives identical results as A*x<b version", {
  ############ indicator function in R: x1 < x2 < x3
  opt <- 3
  k <- c(4, 1, 53)
  mod <- function(x) x[1] < x[2] & x[2] < 1 - x[1] - x[2]
  set.seed(12345)
  expect_silent(mcmc_r <- sampling_nonlinear(k, opt, mod,
    M = 2000,
    progress = FALSE
  ))
  expect_true(all(apply(mcmc_r, 1, mod)))
  expect_silent(bf_r <- bf_nonlinear(k, opt, mod, M = 10000, progress = FALSE))

  ############# A*x < b
  A <- matrix(c(
    1, -1,
    1, 2
  ), 2, 2, TRUE)
  b <- c(0, 1)
  mcmc_Ab <- sampling_multinom(k, opt, A, b, M = 10000, progress = FALSE)
  expect_silent(bf_Ab <- bf_multinom(k, opt, A, b, M = 100000, progress = FALSE))

  expect_gt(ks.test(mcmc_Ab[, 1], mcmc_r[, 1])$p, .005)
  expect_gt(ks.test(mcmc_Ab[, 2], mcmc_r[, 2])$p, .005)

  expect_equal(bf_Ab[, 1], bf_r[, 1], tolerance = max(bf_r[, 2]) * 3)


  ############ indicator function in C++: x1 < x2 < x3
  skip_on_cran()
  skip_if_not_installed("RcppXPtrUtils")

  mod_cpp <- "
    SEXP inside(NumericVector x){
      return wrap(x[0] < x[1] & x[1] < 1-sum(x));
    }" # C++ indices start at 0!
  mod_ptr <- RcppXPtrUtils::cppXPtr(mod_cpp)
  expect_silent(mcmc_cpp <- sampling_nonlinear(k, opt, mod_ptr,
    M = 10000,
    progress = FALSE
  ))
  expect_gt(ks.test(mcmc_Ab[, 1], mcmc_cpp[, 1])$p, .005)
  expect_gt(ks.test(mcmc_Ab[, 2], mcmc_cpp[, 2])$p, .005)

  expect_silent(bf_cpp <- bf_nonlinear(k, opt, mod_ptr, M = 100000, progress = FALSE))
  expect_equal(bf_Ab[, 1], bf_cpp[, 1], tolerance = max(bf_cpp[, 2]) * 3)
})


# test_that("nonlinear Gibbs/BF gives correct results for contigency tables",{
#
# })
