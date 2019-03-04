library("testthat")


# vertex representation (one prediction per row)
V <- matrix(c(
  # strict weak orders
  0, 1, 0, 1, 0, 1,  # a < b < c
  1, 0, 0, 1, 0, 1,  # b < a < c
  0, 1, 0, 1, 1, 0,  # a < c < b
  0, 1, 1, 0, 1, 0,  # c < a < b
  1, 0, 1, 0, 1, 0,  # c < b < a
  1, 0, 1, 0, 0, 1,  # b < c < a

  0, 0, 0, 1, 0, 1,  # a ~ b < c
  0, 1, 0, 0, 1, 0,  # a ~ c < b
  1, 0, 1, 0, 0, 0,  # c ~ b < a
  0, 1, 0, 1, 0, 0,  # a < b ~ c
  1, 0, 0, 0, 0, 1,  # b < a ~ c
  0, 0, 1, 0, 1, 0,  # c < a ~ b

  0, 0, 0, 0, 0, 0   # a ~ b ~ c
), byrow = TRUE, ncol = 6)

# transform to Ab-representation
# Ab <- V_to_Ab(V)

# for CRAN: use pre-computed solution provided by rPorta
Ab <- list(A = structure(c(-1, 0, 0, 0, 0, 0, -1, -1, 0, 0, 0, 1, 0,
                           0, 1, 0, -1, 0, 0, 0, 0, 0, 0, -1, -1, 1, 0, 0, 0, 1, 0, 0, -1,
                           0, 0, 0, 0, 1, -1, 0, 0, -1, 0, 1, 0, 0, 0, 0, -1, 0, 0, -1,
                           0, 0, 1, -1, 0, 0, 1, 0, 0, 0, 0, 0, -1, 0, 0, -1, 1, 0, -1,
                           0, 1, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0, -1, 0, -1, 1, 0, 0),
                         .Dim = c(15L, 6L)),
           b = c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1))

test_that("ML estimation matches for Ab- and V-representation", {

  set.seed(123)
  k <- c(4,1,5,  1,9,0,  7,2,1)
  options <- c(3, 3, 3)

  est_V <- ml_multinom(k = k, options = options, V = V,
                       n.fit = 5, control = list(maxit = 1e6, reltol=.Machine$double.eps^.6),
                       outer.iterations = 1000
                       )


  est_Ab <- ml_multinom(k = k, options = options, A = Ab$A, b = Ab$b,
                        n.fit = 5, control = list(maxit = 1e6, reltol=.Machine$double.eps^.6),
                        outer.iterations = 1000
                        )

  expect_equal(unname(est_V$p), est_Ab$par, tolerance = .0001)
})

test_that("if unconstrained MLE == constrained MLE", {

  set.seed(123)
  n <- 100
  m <- 5
  options = rep(3,m)
  V <- round(rpdirichlet(n, rep(1,3*5), options = options))
  alpha <- c(rpdirichlet(1, rep(1, n), n, drop_fixed = FALSE))
  p <- add_fixed( c(alpha %*% V), options)
  k <- c(round(p*10000))
  n <- c(tapply(k, rep(1:length(options), options), sum))
  mle_unconstr <- drop_fixed(k / n, options)

  t <- system.time(
    est_V <- ml_multinom(k = k, options = options, V = unique(V),
                         n.fit = 5, control = list(maxit = 1e7, reltol=.Machine$double.eps^.5),
                         outer.iterations = 1000, progress = FALSE))["elapsed"]
  expect_lt(t, 1)  # should be fast!
  expect_equal(est_V$p, mle_unconstr, tolerance = .0001)
  expect_equal(est_V$value, multinomineq:::loglik_multinom(mle_unconstr, k, options))

})
