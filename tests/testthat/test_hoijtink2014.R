library("testthat")
library("multinomineq")


# Model: p1 > p2 > ... > p6
k <- c(1, 1, 1, 0, 0, 0)
A <- -diag(6)[-6, ] + cbind(0, diag(6)[-6, -6])
b <- rep(0, 5)
cbind(A = A, b = b)


test_that("results match with cognitive diagnostic testing (Hoijtink 2014)", {
  set.seed(123)
  cnt <- count_binom(k, 1, A, b, M = 5e5)
  expect_equal(attr(cnt, "proportion"), .011, tolerance = .001)

  # .011 * factorial(6) / ((1-.011) * (1-1/factorial(6)))

  bf <- count_to_bf(cnt, exact_prior = 1 / factorial(6))
  expect_equal(bf[3, 1], 7.76, tolerance = 3 * bf[3, 2]) # cf. p. 32, table 6


  # function for log(1-exp(x))
  xx <- log(runif(10))
  expect_equal(multinomineq:::log1mexp(xx), log(1 - exp(xx)))
})
