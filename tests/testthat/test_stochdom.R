
# test_that("stochastic dominance is tested correctly", {
#
#   bfs <- stochdom_bf(x1, x2, order = "<", cmin=1000, M=5000, steps = 1)$bf
#   bfs
#   j <- rep(1/nrow(b12$k), 2*(nrow(b12$k)))
#   s <- ml_multinom(c(b12$k), rep(nrow(b12$k), 2), b12$A, b12$b, n.fit = 1)
#   post <- count_multinom(c(b12$k), rep(nrow(b12$k), 2), b12$A, b12$b, M=500, prior=j, cmin=10)
#   post
#   count_multinom(c(b12$k), rep(nrow(b12$k), 2), b12$A, b12$b, M=5e5, prior=j)
#
#   prior <- count_multinom(0, rep(nrow(b12$k), 2), b12$A, b12$b,M=1000, prior=j,  cmin=10)
#   prior
#   count_multinom(0, rep(nrow(b12$k), 2), b12$A, b12$b,M=2e5, prior=j)
#   count_to_bf(post,prior)
#
# })

