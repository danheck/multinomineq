library(testthat)
set.seed(123)

M <- 2e5
p <- seq(0,1,.05)
quantile_ss <- function(x, y, p = seq(0,1,.05)){
  sum((quantile(x, p) - quantile(y, p))^2)
}
sel_Ab <- function(X, A, b)
  apply(X, 1, function(x) all(A %*% x <= b))
test_Ab <- function(X, A, b){
  all(sel_Ab(X,A,b))
}
test_sum1 <- function(X, options){
  all.equal(c(apply(X, 1, function(x)
    tapply(x, rep(1:length(options), options), sum))),
    rep(1, length(options)*nrow(X)))
}
test_sum1_free <- function(X, options){
  all(apply(X, 1, function(x)
    tapply(x, rep(1:length(options), options-1), sum)) <= 1)
}

test_that("truncated beta is correct", {
  bmin <- .236
  bmax <- 0.779
  a <- 14
  b <- 4
  x <- replicate(5000, multinomineq:::rbeta_trunc(a, b, bmin, bmax))
  y <- rbeta(10000, a, b)
  x2 <- y[y > bmin & y < bmax]
  expect_equal(quantile(x, p), quantile(x2, p), tol = .02)
  # qqplot(x, x2)
})

test_that("Gibbs sampling gives same results as accept-reject for binomial", {

  ################### prior: uniform
  k <- n <- rep(0, 5)
  A <- Ab_multinom(rep(2, length(n)))$A
  b <- Ab_multinom(rep(2, length(n)))$b
  X <- sampling_binom(k, n, A, b, M = M)
  expect_true(test_Ab(X, A, b))
  expect_equal(colMeans(X), 1/rep(2, ncol(A)), tol = .02)

  ################### posterior constrained (vs. accept-reject)
  k <- c(15,9,5,2)
  n <- c(20, 25, 15, 20)
  A <- matrix(c(1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 2, 1, 0,
                0, -1, 0,0,
                0, 0, 1 ,0,
                0,0,  1,1,
                0,0,0,  2),
              ncol = length(k), byrow = TRUE)
  b <- c(.7, .25, .25, 1, -.2, .8, 1,1)
  start = c(0.1048014, 0.2141486, 0.1246359, 0.3551680)

  X1 <- sampling_binom(k, n, A, b, M = M, start = start)
  expect_true(test_Ab(X1, A, b))
  # accept/reject
  X <- matrix(rbeta(10*M*length(k), k+1,n-k+1), ncol = length(k), byrow=TRUE)
  sel <- apply(X, 1, function(x) all(A %*% x <= b))

  # par(mfrow=c(1,2))
  # for (i in 1:ncol(A)){
  #   plot(density(X1[,i]))
  #   lines(density(X[sel,i]), col = 2)
  #   qqplot(X1[,i], X[sel,i])
  # }
  expect_equal(colMeans(X1), colMeans(X[sel,]), tol = .01)
  expect_equal(unname(apply(X1, 2, sd)), apply(X[sel,], 2, sd), tol = .01)
  for (i in 1:ncol(A))
    expect_equal(quantile(X1[,i], p), quantile(X[sel,i], p), tol = .04)
})


test_that("Gibbs sampling for multinomial [prior: k=0]", {

  options <- 4
  A <- matrix(1, 1, options - 1)
  b <- 1
  k <- rep(0, sum(options))
  X <- sampling_multinom(k, options, A, b, M = M)
  X2 <- multinomineq:::rdirichlet(M, k + 1)
  expect_true(test_sum1_free(X[1:1000,], options))
  expect_equal(unname(colMeans(X)), 1/rep(options, options - 1), tol = .005)
  for (i in 1:ncol(X))
    expect_equal(quantile_ss(X[,i],X2[,i],p), 0, tol = .1)

  # binary and ternary choice:
  options <- c(2,     3,      4)

  ################### unconstrained sampling from prior:
  A <- Ab_multinom(options)$A
  b <- Ab_multinom(options)$b
  k <- rep(0, sum(options))
  X <- sampling_multinom(k, options, A, b, M = M)
  expect_true(test_sum1_free(X[1:1000,], options))
  expect_equal(unname(colMeans(X)), 1/rep(options, options - 1), tol = .005)
  expect_equal(unname(quantile(X[,1], p)), p, tol = .05)

  ################### uniform prior comparison to accept/reject
  # columns:   (a1,  b1,b2)
  A <- matrix(c(1, 0, 0, 0,0,0,  # a1 < .20
                0, 2, 1, 0,0,0,  # 2*b1+b2 < 1
                0, -1, 0,0,0,0,  # b1 > .2
                0, 0, 1 ,0,0,0,
                0,0,0,  1,1,0,
                0,0,0, 0, 2,3),  # b2 < .4
              ncol = sum(options-1), byrow = TRUE)
  b <- c(.2, 1, -.2, .4, 1,1)
  k <- rep(0, sum(options))
  start <- c(0.1871862, 0.2881284, 0.3919744, 0.4760341, 0.2890125, 0.1379431)
  X1 <- sampling_multinom(k, options, A, b = b, M = 2000, start = start)
  expect_true(test_sum1_free(X1, options))
  expect_true(test_Ab(X1, A, b))
  X2 <- rpdirichlet(M*5, k+1, options, drop_fixed = TRUE)
  sel <- sel_Ab(X2, A, b)
  for (i in 1:ncol(X1)){
    # qqplot(X1[,i],X2[sel,i])
    # abline(0,1,col=2)
    expect_equal(quantile_ss(X1[,i],X2[sel,i],p), 0, tol = .05)
  }
})

test_that("Gibbs sampling for multinomial [posterior]", {

  ################### one constrained DIrichlet / multinomial
  options <- c(4)
  A <- matrix(c(1,1,1,
                1,3,0,
                -1, 0, 0,
                0, -1, 0,
                0, 3, -1),
              ncol = sum(options-1), byrow = TRUE)
  b <- c(.7, 1.5, -.2, -.3, 1)
  k <- c(8, 3, 0, 12)
  n <- multinomineq:::sum_options(k, options)
  expect_equal(c(n), unname(rep(tapply(k, rep(1:length(options), options), sum), options)))
  start <- c(0.2531039, 0.3094770, 0.0494943)
  X <- sampling_multinom(k, options, A, b, M = M, start=start)
  expect_true(test_sum1_free(X[1:100,], options))
  expect_true(test_Ab(X[1:1000,], A, b))
  # accept-reject
  X1 <- rpdirichlet(M*5, k + 1, options, drop_fixed = TRUE)
  sel <- sel_Ab(X1, A, b)
  # par(mfrow=c(2,2))
  # for (i in 1:ncol(A)){
  #   plot(density(X[,i]), xlim = 0:1, main = i)
  #   lines(density(X1[sel,i]), col = 2)
  #   qqplot(X[,i], X1[sel,i])
  #   abline(0,1,col=2)
  # }
  expect_equal(unname(colMeans(X)), colMeans(X1[sel,]), tol = .01)
  expect_equal(unname(apply(X, 2, sd)), apply(X1[sel,], 2, sd), tol = .005)
  for (i in 1:ncol(A))
    expect_equal(quantile_ss(X[,i],X1[sel,i], p), 0, tol = .001)
})

test_that("Gibbs sampling for product-multinomial [posterior]", {

  ################### posterior   product-multinomial/dirichlet
  options <- c(2,     3,      4)
  A <- matrix(c(1, 0, 0, 0,0,0,  # a1 < .70

                0, -1, 0,0,0,0,  # b1 > .2
                0, 0, 1 ,0,0,0,
                0, 0, -1, 0,0,0,

                0,0,0,  -1,0,0,
                0,0,0,  0,1,1,
                0,0,0,  1,0,1,

                1,1,1,  1,1,1),  # b2 < .4
              ncol = sum(options-1), byrow = TRUE)
  b <- c(.5, -.2, .8, -.1, -.2, 1, .53, 1.8)
  start <- c(0.1871862, 0.2881284, 0.3919744, 0.3760341, 0.2890125, 0.1379431)
  k       <- c(15,9,   5,2,17, 5,1,8,20)
  X <- sampling_multinom(k, options, A, b, M = M, start = -1)
  expect_true(test_sum1_free(X[1:100,], options))
  expect_true(test_Ab(X[1:100,], A, b))
  # accept-reject
  X1 <- rpdirichlet(M*2, k + 1, options, drop_fixed = TRUE)
  sel <- sel_Ab(X1, A, b)
  cnt <- count_multinom(k, options, A, b, M = 5e5)
  expect_equal(mean(sel), attr(cnt, "proportion"), tol = 5 * attr(cnt, "se")) ## integral
  # par(mfrow=c(2,2))
  # for (i in 1:ncol(A)){
  #   plot(density(X[,i]), xlim = 0:1, main = i)
  #   lines(density(X1[sel,i]), col = 2)
  #   qqplot(X[,i], X1[sel,i])
  # }
  expect_equal(unname(colMeans(X)), colMeans(X1[sel,]), tol = .005)
  expect_equal(unname(apply(X, 2, sd)), apply(X1[sel,], 2, sd), tol = .005)
  for (i in 1:ncol(A))
    expect_equal(quantile_ss(X[,i],X1[sel,i], p), 0, tol = .03)

  rm(X); rm(X1)
  gc()

  ############## swop polytope
  data(swop5)
  options <- rep(3, 10)
  k <- rep(0, 30)
  X <- rpdirichlet(M, k + 1, options, drop_fixed = TRUE)
  int <- multinomineq:::count_samples(X, swop5$A, swop5$b)/(M)
  cnt <- count_multinom(k, options, swop5$A, swop5$b, M=M)
  expect_equal(int, attr(cnt, "proportion"), tol = 5*attr(cnt, "se"))

  # X2 <- sampling_multinom(k, swop5$options, swop5$A, swop5$b,
  #                         M = 500, start = swop5$start)
  # sel <- apply(X, 1, inside, A =swop5$A, b=swop5$b)
  # rbind(colMeans(X[sel,]), colMeans(X2))
})
