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

set.seed(124356)

test_that("counting methods work", {
  cnt <- stratsel:::count_samples(X, A, b)
  expect_equal(cnt, cntr(X, A, b))
  expect_equal(cnt / M, V, tol = .005)

  c2 <- count_binom(0, 0, A=A, b=b, M = M, batch = M/10)
  expect_equal(attr(c2, "integral"), V, tol = .005)
  c3 <- count_binom(0, 0, A=A, b=b, M = M, batch = M)
  expect_equal(attr(c3, "integral"), V, tol = .005)

  c4 <- count_binom(0, 0, A=A, b=b, M = 1e5, batch = 10000, steps = 2:3)
  c4
  expect_equal(attr(c4, "integral"), V, tol = .0005)
})

test_that("counting is equivalent for: A/b-method and V-method", {
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

  xin <- c(.1, .2, .45)  # inside
  expect_true(inside(xin, A, b))
  expect_true(inside(xin, V = V))
  xout <- c(.4, .1, .55)  # outside
  expect_false(inside(xout, A, b))
  expect_false(inside(xout, V = V))
  expect_equal(inside(rbind(xin, xout), A, b), c(TRUE, FALSE))
  expect_equal(inside(rbind(xin, xout), V = V), c(TRUE, FALSE))


  V <- matrix(c(1,0,0,0,0,1,
                0,1,1,0,0,0,
                1,0,0,1,0,0,
                1,1,1,1,1,0,
                1,0,0,1,1,0,
                0,1,1,1,0,1,
                0,0,1,1,1,1,
                0,0,0,0,0,0,
                0,1,1,1,0,1,
                1,1,1,1,1,1,
                0,1,0,1,0,1,
                0,0,0,0,0,1), ncol=6, byrow=TRUE)
  tmp <- V_to_Ab(V)
  A <- tmp$A
  b <- tmp$b

  X <- matrix(runif(2000 * ncol(V)), ncol = ncol(V))
  vertex <- stratsel:::inside_V(X, V)
  ineq <- apply(X, 1, function(x) all(A %*% x <= b))
  expect_equal(vertex, ineq)
  table(data.frame(vertex, ineq))

})

test_that("transformation of A/b to V representation works",{
  skip_on_cran()
  skip_if_not_installed("rPorta")


})
