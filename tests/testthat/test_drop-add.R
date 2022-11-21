library("testthat")


test_that("fixed probabilities/numbers are dropped correctly with 'drop_fixed'", {
  # vector
  k <- c(1, 4, 2, 3, 1, 4)
  expect_silent(k1 <- drop_fixed(k, rep(2, 3)))
  expect_equal(unname(k1), c(1, 2, 1))
  expect_silent(k2 <- drop_fixed(k, 2))
  expect_equal(k1, k2)

  expect_error(drop_fixed(c(1, 3, 2, 4), c(2, 3)))
  expect_error(drop_fixed(c(1, 3, 2, 4), c(2, 3, 2)))
  expect_error(drop_fixed(c(1, 3, 2, 4), 5))

  # matrix
  k <- matrix(c(1, 4, 2, 3, 1, 4), 1)
  expect_silent(k1 <- drop_fixed(k, rep(2, 3)))
  expect_equal(unname(k1), matrix(c(1, 2, 1), 1))
  expect_silent(k2 <- drop_fixed(k, 2))
  expect_equal(k1, k2)

  expect_silent(k3 <- drop_fixed(k, c(2, 4)))
  expect_equal(unname(k3), k[, c(1, 3:5), drop = FALSE])

  expect_error(drop_fixed(k, c(2, 3)))
  expect_error(drop_fixed(k, 4))
})
