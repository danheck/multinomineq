library(multinomineq)
library(testthat)
# sampling vertex polytope
options <- c(3,2)
V <- matrix(c( 0, 0, 0,
               0, 0,.5,
               0,.5,.5,
               .5,.5,.5), ncol = 3, byrow = TRUE)
# V <- matrix(c( 0, 0, 0,
#                0, 0,1,
#                0,1,0,
#                1,0,0,
#                1,0,1,
#                0,1,1), ncol = 3, byrow = TRUE)
V <- t(apply(V, 1, multinomineq:::add_fixed, options = options))
k <- c(3,4,2,  1,5)
Ab <- V_to_Ab(V[,c(1:2,4)])

test_that("posterior sampling for vertex representation works", {

  ###### check: accept-reject
  u <- rpdirichlet(1e5, 1 + k, options)
  sel <- apply(u[,c(1:2,4)], 1, inside, A=Ab$A, b = Ab$b)
  s.AR <- u[sel,c(1:2,4),drop = FALSE]
  s.V <- sampling_V(k, options, V, M = 5000, start = s.AR[nrow(s.AR),])

  # check whether all samples are inside polytope
  expect_true(all(apply(s.V, 1, inside, A = Ab$A, b = Ab$b)))

  for(i in 1:3){
    expect_gte(suppressWarnings(ks.test(s.V[,i], s.AR[,i])$p), .01)
    qqplot(s.V[,i], s.AR[,i])
    abline(0, 1, col=2)
  }

})



