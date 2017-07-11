d <- 4
A <- matrix(c(1, -1, 0, 0,
              0, 1, -1, 0,
              0, 0, 1, -1,
              0, 0, 0, 1), ncol = d, byrow = TRUE)
b <- c(0,0,0,.5)
V <- .5^d / factorial(d)

M <- 50000
X <- matrix(runif(M*d), M)

cntr <- function(X, A, b)
  sum(apply(X, 1, function(x) all(A %*% x <= b)))

test_that("counting methods work", {
  cnt <- stratsel:::count_samples(X, A, b)
  expect_equal(cnt, cntr(X, A, b))
  expect_equal(cnt / M, V, tol = .002)

  c2 <- count_polytope(A, b, M = M, batch = M/10)
  expect_equal(c2$integral, V, tol = .002)
  c3 <- count_polytope(A, b, M = M, batch = M)
  expect_equal(c3$integral, V, tol = .002)

})

