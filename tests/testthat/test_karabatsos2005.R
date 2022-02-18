library(testthat)

# compare results to
# Karabatsos, G. (2005). The exchangeable multinomial model as an approach to testing deterministic axioms of choice and measurement. Journal of Mathematical Psychology, 49(1), 51-69. https://doi.org/10.1016/j.jmp.2004.11.001



test_that("results match to those of Karabatsos (2005)", {
  # example from appendix:
  H0 <- list(pattern =  1, c = .5, ordered = FALSE, prior = c(.5,.5))  # predict: non-violations (k = n)
  H1 <- list(pattern = -1, c = .5, ordered = FALSE, prior = c(.5,.5))
  expect_silent(bf <- strategy_postprob(cbind(0:5), cbind(rep(10,6)), list(H0, H1)))
  expect_equal(bf[,1], c(.00016, .00369, .02604, .10202, .26483,.5), tol = .00001)

  # Table 4: Gambles 2-4 with reference prior (1/2, 1/2):
  H0 <- list(pattern = c(1), c = .5, ordered = FALSE, prior = c(.5,.5))
  H1 <- list(pattern = c(-1),c = .5, ordered = FALSE, prior = c(.5,.5))
  expect_silent(bf <- strategy_postprob(cbind(c(15, 17, 18)), cbind(rep(31,3)), list(H0, H1)))
  expect_equal(bf[,1] / bf[,2], c(.75, 2.4, 4.4), tol = .1)

  # Gambles 2-4 with different prior:  [Karabatsos: see appendix for syntax]
  s <- c(15, 17, 18)
  prior <- pbeta(1, 4, 1) - pbeta(.5, 4, 1)
  post  <- pbeta(1, 4+s, 1+31-s) - pbeta(.5, 4+s, 1+31-s)
  bf_0u <- post/prior
  bf_1u <- (1-post)/(1-prior)

  H0 <- list(pattern = c(1), c = .5, ordered = FALSE, prior = c(1, 4)) # tau_violations, tau_ok
  en <- list(pattern = c(1), c =  1, ordered = FALSE, prior = c(1, 4))
  expect_silent(bf <- strategy_postprob(cbind(c(15, 17, 18)), cbind(rep(31,3)), list(H0, en)))
  expect_equal(bf[,1] / bf[,2], bf_0u)

  H0 <- list(pattern = c(1), c = .5, ordered = FALSE, prior = c(1, 4)) # tau_violations, tau_ok
  H1 <- list(pattern = c(-1),c = .5, ordered = FALSE, prior = c(4, 1))
  expect_silent(bf <- strategy_postprob(cbind(c(15, 17, 18)), cbind(rep(31,3)), list(H0, H1)))

  expect_equal(bf[,1] / bf[,2], bf_0u/bf_1u)



  # assumption of independence for gambles 1-5
  H0 <- list(pattern = 1:5,   c = .5, ordered = FALSE, prior = c(.5,.5))
  H1 <- list(pattern = -(1:5),c = .5, ordered = FALSE, prior = c(.5,.5))
  expect_silent(bf <- strategy_postprob(c(23, 15, 17, 18, 21), rep(31,5),
                                list(H0, H1)))
  expect_equal(bf[1] / bf[2], 109181, tol = 1)

})
