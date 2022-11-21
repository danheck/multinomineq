library("testthat")

# simulated data in Table 1
klaassen2018 <- matrix(c(
  7, 5, 4, 1,
  7, 2, 5, 1,
  3, 4, 6, 1
), 3, byrow = TRUE)
A1 <- matrix(c(
  -1, 1, 0, 0,
  0, -1, 1, 0,
  0, 0, -1, 1
), 3, byrow = TRUE)
b1 <- rep(0, 3)

A2 <- matrix(c(-1, -1, 1, 1), 1)
b2 <- 0

set.seed(124342)

test_that("Klaassen Table3: Hyp. 1 results match", {
  expect_silent(cb1 <- lapply(data.frame(t(klaassen2018)), count_binom,
    n = 7, A = A1, b = b1, cmin = 1e4, M = 1e4, progress = FALSE
  ))
  bf_1u <- lapply(cb1, count_to_bf, exact_prior = 1 / factorial(4))

  bfs <- sapply(bf_1u, "[", "bf_0u", c("bf", "se"))
  expect_equal(unname(bfs["bf", ]), c(13.16, 1.40, .24), tol = max(bfs["se", ]) * 3)
})

test_that("Klaassen Table3: Hyp. 2 results match", {
  expect_silent(cb2 <- lapply(data.frame(t(klaassen2018)), count_binom,
    n = 7, A = A2, b = b2, M = 3e5, progress = FALSE
  ))
  bf_2u <- lapply(cb2, count_to_bf, exact_prior = .5)

  bfs <- sapply(bf_2u, "[", "bf_0u", c("bf", "se"))
  expect_equal(unname(bfs["bf", ]), c(2.00, 1.79, 1.01), tol = max(bfs["se", ] * 3, .03))
})
