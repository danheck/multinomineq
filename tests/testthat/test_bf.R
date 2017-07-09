# models of Hilbig & Moshagen (2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
preds <- lapply(list(TTB=TTB, WADD=WADD, WADDprob=WADDprob,
                     EQW=EQW, GUESS=GUESS, baseline=baseline),
                as_strategy)
preds$baseline$c <- .5
preds$baseline$ordered <- FALSE


k <- c(0, 0, 0)
k_false <- c(99,1,1)
n <- rep(32,3)
prior <- c(1, 1)

test_that('Bayes factor works as expected', {
  # guessing
  expect_equal(compute_marginal(k, n, as_strategy(GUESS)),
               sum(lchoose(n, k)) + sum(n)*log(.5))

  # TTB
  s1 <- sum(k) + prior[1]
  s2 <- sum(n-k) + prior[2]
  expect_equal(compute_marginal(k, n, as_strategy(TTB)),
               sum(lchoose(n, k)) +
                 pbeta(.5, s1, s2, log = TRUE) +
                 lbeta(s1, s2) - log(.5))

  # TTB is best model:
  expect_named(margs <- select_bf(k, n, preds), names(preds))
  expect_gt(margs[1], max(margs[-c(1)]))

  # WADD is best model:
  expect_named(margs <- select_bf(c(0,n[3],0), n, preds), names(preds))
  expect_gt(margs[2], max(margs[-c(2)]))

  # errors
  expect_error(compute_marginal(k_false, n, as_strategy(TTB)))
})
