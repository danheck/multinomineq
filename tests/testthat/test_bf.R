# models of Hilbig & Moshagen (2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
attr(WADDprob, "ordered") <- TRUE
preds <- list(TTB, WADD, WADDprob, EQW, GUESS, baseline)

b <- c(0, 0, 0)
b_false <- c(99,1,1)
n <- rep(32,3)
prior <- c(1, 1)

test_that('Bayes factor works as expected', {
  # guessing
  expect_equal(compute_marginal(b, n, GUESS),
               sum(lchoose(n, b)) + sum(n)*log(.5))

  # TTB
  s1 <- sum(b) + prior[1]
  s2 <- sum(n-b) + prior[2]
  expect_equal(compute_marginal(b, n, TTB),
               sum(lchoose(n, b)) +
                 pbeta(.5, s1, s2, log = TRUE) +
                 lbeta(s1, s2) - log(.5))

  # TTB is best model:
  expect_silent(margs <- select_bf(b, n, preds))
  expect_gt(margs[1], max(margs[-c(1)]))

  # WADD is best model:
  expect_silent(margs <- select_bf(c(0,n[3],0), n, preds))
  expect_gt(margs[2], max(margs[-c(2)]))

  # errors
  expect_warning(compute_marginal(b_false, n, TTB))
})
