library(multinomineq)
library(testthat)
# sampling vertex polytope
options <- c(3,2)
V <- matrix(c( 0, 0, 0,
               0, 0,.5,
               0,.5,.5,
               .5,.5,.5), ncol = 3, byrow = TRUE)
k <- c(3,4,2,  1,5)
Ab <- list(A = structure(c(-1, 0, 1, 0, 0, 1, -1, 0, 0, -1, 0, 2), .Dim = 4:3),
           b = c(0, 0, 0, 1))

test_that("posterior sampling for vertex representation works", {

  set.seed(123)

  ###### check: accept-reject
  u <- rpdirichlet(1e5, 1 + k, options, drop_fixed = TRUE)
  sel <- apply(u, 1, inside, A=Ab$A, b = Ab$b)
  s.AR <- u[sel,,drop = FALSE]
  s.V <- sampling_multinom(k, options, V = V, M = 2000, start = s.AR[nrow(s.AR),])

  # check whether all samples are inside polytope
  expect_true(all(apply(s.V, 1, inside, A = Ab$A, b = Ab$b)))

  for(i in 1:3){
    expect_gte(suppressWarnings(ks.test(s.V[,i], s.AR[,i])$p), .01)
    qqplot(s.V[,i], s.AR[,i])
    abline(0, 1, col=2)
  }

})



