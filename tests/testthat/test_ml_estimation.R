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
