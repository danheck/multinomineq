library(testthat)
test_that("transformation of A/b to V representation works",{
  skip_on_cran()
  skip_if_not_installed("rPorta")

  A <- matrix(c(1,-1, 0,
                0, 1,-1,
                0, 0, 1), ncol = 3, byrow = TRUE)
  b <- c(0, 0, .50)
  V <- matrix(c( 0, 0, 0,
                 0, 0,.5,
                 0,.5,.5,
                 .5,.5,.5), ncol = 3, byrow = TRUE)

  # tmp <- V_to_Ab(V)

  expect_equal(Ab_to_V(A, b), V)
})
