
dim <- 4
k <- c(0, 5, 9, 12)
n <- rep(30, 4)
# linear order constraint:
# theta1 < theta2 < .... < .50
A <- matrix(c(1, -1, 0, 0,
              0, 1, -1, 0,
              0, 0, 1, -1,
              0, 0, 0, 1),
            ncol = 4, byrow = TRUE)
b <- c(0, 0, 0, .50)

set.seed(123)
test_that("encompassing Bayes factors returns correct results", {

  # exact volume of linear order constraint
  volume <- 1/factorial(dim)/2^dim

  # standard encompassing BF
  expect_silent(bf1 <- compute_bf(k, n, A, b, c(1, 1), M = 1e5))
  expect_named(bf1, c("bf", "se", "posterior", "prior", "M"))
  expect_named(bf1$bf, c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"))
  expect_named(bf1$se, c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"))
  expect_equal(5.3, unname(bf1$bf["log_bf_0e"]), tolerance = .5)
  expect_equal(log(volume), log(bf1$prior/bf1$M), tolerance = .3)

  # invariance under reordering of order inequalities
  expect_silent(res <- sort_inequalities(A, b, M = 2e5))
  A2 <- res$A
  b2 <- res$b
  bf3 <- compute_bf(k, n, A2, b2, c(1, 1), M = 1e5)
  bf4 <- compute_bf(k, n, A2, b2, c(1, 1), M = c(2e5, 2e4, 2e4), steps = c(1,3))
  expect_equal(5.3, unname(bf3$bf["log_bf_0e"]), tolerance = .2)
  expect_equal(5.3, unname(bf4$bf["log_bf_0e"]), tolerance = .2)
})

