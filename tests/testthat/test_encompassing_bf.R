
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

  set.seed(1234)

  # exact volume of linear order constraint
  volume <- 1/factorial(dim)/2^dim
  vol_est <- count_binom(k = 0, n = 0, A = A , b = b, M = 2e5)
  expect_equal(volume, attr(vol_est, "proportion"), tolerance = 10 * attr(vol_est, "se"))

  # posterior
  post <- count_binom(k, n, A, b, prior = c(1, 1), M = 2e5)
  bf <- count_to_bf(post, vol_est, log = TRUE)
  expect_equal(5.3, bf["bf_0u","bf"], tolerance = 10 * bf["bf_0u","se"])
  bf_exact <- count_to_bf(post, exact_prior = volume, log = TRUE)
  expect_equal(5.3, bf_exact["bf_0u", "bf"], tolerance = 10 *  bf_exact["bf_0u", "se"])

  # standard encompassing BF
  expect_silent(bf1 <- bf_binom(k, n, A, b, prior = c(1, 1), M = 2e5, log = TRUE))
  expect_equal(5.3, bf1["bf_0u","bf"], tolerance = 10 * bf1["bf_0u","se"])#


  # invariance under reordering of order inequalities
  expect_silent(res <- sort_Ab(A, b, M = 2e5))
  A2 <- res$A
  b2 <- res$b
  bf3 <- bf_binom(k, n, A2, b2, prior = c(1, 1), M = 1e5, log = TRUE)
  bf4 <- bf_binom(k, n, A2, b2, prior = c(1, 1), M = c(2e5, 1e4, 1e4), steps = c(1,3), log = TRUE)
  expect_equal(5.3, bf3["bf_0u","bf"], tolerance = .2)
  expect_equal(5.3, bf4["bf_0u","bf"], tolerance = .2)
})

