
dim <- 4
k <- c(0, 5, 9, 12)
n <- rep(30, dim)
# linear order constraint:
# x1 < x2 < x3 < x4 < .50
A <- matrix(c(1, -1, 0, 0,
              0, 1, -1, 0,
              0, 0, 1, -1,
              0, 0, 0, 1),
            ncol = dim, byrow = TRUE)
b <- c(0, 0, 0, .50)

# set.seed(123)
test_that("encompassing Bayes factors returns correct results", {

  # exact volume of linear order constraint
  volume <- 1/factorial(dim)/2^dim

  # standard encompassing BF
  expect_silent(bf1 <- bf_binom(k, n, A, b, c(1, 1)))
  # expect_named(bf1, list(c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"), c("BF", "SE")))
  expect_equal(5.3, bf1["log_bf_0e","BF"], tolerance = .5)

  # invariance under reordering of order inequalities
  expect_silent(res <- sort_Ab(A, b, M = 2e5))
  A2 <- res$A
  b2 <- res$b
  bf3 <- bf_binom(k, n, A2, b2, c(1, 1), M = 1e5)
  bf4 <- bf_binom(k, n, A2, b2, c(1, 1), M = c(2e5, 1e4, 1e4), steps = c(1,3))
  expect_equal(5.3, bf3["log_bf_0e","BF"], tolerance = .2)
  expect_equal(5.3, bf4["log_bf_0e","BF"], tolerance = .2)
})

