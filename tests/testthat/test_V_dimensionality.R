

### ternary choice (Regenwetter & Davis-Stober, 2012)
# choice options:  {prefer_a, indifferent, prefer_b}
library("testthat")
# column order:    (a>b,b>a,  a>c,c>a,  b>c,c>b)
# with:            i>j = 1  <=> utility(i) > utility(j)
V <- matrix(c(
  # strict weak orders
  0,1,0, 0,1,0, 0,1,0,  # a < b < c
  1,0,0, 0,1,0, 0,1,0,  # b < a < c
  0,1,0, 0,1,0, 1,0,0,  # a < c < b
  0,1,0, 1,0,0, 1,0,0,  # c < a < b
  1,0,0, 1,0,0, 1,0,0,  # c < b < a
  1,0,0, 1,0,0, 0,1,0,  # b < c < a

  0,0,1, 0,1,0, 0,1,0,  # a ~ b < c
  0,1,0, 0,0,1, 1,0,0,  # a ~ c < b
  1,0,0, 1,0,0, 0,0,1,  # c ~ b < a
  0,1,0, 0,1,0, 0,0,1,  # a < b ~ c
  1,0,0, 0,0,1, 0,1,0,  # b < a ~ c
  0,0,1, 1,0,0, 1,0,0,  # c < a ~ b

  0,0,1, 0,0,1, 0,0,1   # a ~ b ~ c
), byrow = TRUE, ncol = 9)
options <- rep(3,3)

test_that("dimensionality of V works", {
  expect_silent(V_free <- drop_fixed(V, options))
  expect_is(V_free, "matrix")
  expect_equal(dim(V_free), c(nrow(V), sum(options - 1)))

  # only with Porta:
  # V_to_Ab(V_free)

  expect_silent(p <- find_inside(V=V, random = TRUE))
  expect_true(inside(p, V = V))
  expect_true(inside(drop_fixed(p, options), V=V_free))
  expect_silent(p_free <- find_inside(V=V_free, random = TRUE))
  expect_true(inside(p_free, V = V_free))
  expect_true(inside(add_fixed(p_free, options), V=V))


  # undebug(multinomineq:::sampling_V)
  expect_silent(pp <- sampling_multinom(k = c(4,2,3,  19,4,2,  2,15,10),
                                        options = rep(3,3), V = V_free,
                                        M = 100, progress = FALSE))
  expect_true(all(inside(pp, V = V_free)))
  expect_true(all(inside(add_fixed(pp, options), V = V)))

  ################### here: V with fixed dimensions!
  # library(Rglpk)
  # d <- c(1,0,0, 0,0,0, 0,0,0)
  # # variables: c(lambda, alpha_1, ..., alpha_M)
  # obj <- c(1, rep(0, nrow(V)))
  #
  # mat <- rbind(rep(c(0, 1), c(1, nrow(V))), # sum(alpha) = 1
  #              cbind(d, - t(V)))            # intersection p+lambda*d = sum(alpha*V)
  # dir <- c(rep("==", 1 + ncol(V)))
  # rhs <- c(1, - p)
  # lp <- Rglpk_solve_LP(obj, mat, dir, rhs, max = TRUE)
  # target <- p + (lp$optimum+.000000001) * d
  # inside(target, V = V)
  # p[d] + lp$optimum
  # target
  # p

  ####  alpha >= 0 by default in Rglpk!
  # M <- nrow(V)
  # alpha.sum <- rep(c(0, 1), c(1, M))
  # line.intersect <- cbind(d, - t(V))
  # alpha.pos <- cbind(0, diag(M))
  # dir <- c(rep(">=", M), rep("==", 1 + ncol(V)))
  # rhs <- c(rep(0, M), 1, p)
  # mat <- rbind(alpha.pos, alpha.sum, intersect)
})
