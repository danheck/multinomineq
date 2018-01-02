library("testthat")
# models of Hilbig & Moshagen (2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
preds <- lapply(list(TTB=TTB, WADD=WADD, WADDprob=WADDprob,
                     EQW=EQW, GUESS=GUESS, baseline=baseline),
                multinomineq:::as_strategy)
preds$baseline$c <- .5
preds$baseline$ordered <- FALSE


k <- c(0, 0, 0)
k_false <- c(99,1,1)
n <- rep(32,3)
prior <- c(1, 1)

test_that('Bayes factor works as expected', {
  # guessing
  expect_equal(multinomineq:::strategy_marginal(k, n, multinomineq:::as_strategy(GUESS)),
               sum(lchoose(n, k)) + sum(n)*log(.5))

  # TTB
  s1 <- sum(k) + prior[1]
  s2 <- sum(n-k) + prior[2]
  expect_equal(strategy_marginal(k, n, multinomineq:::as_strategy(TTB)),
               sum(lchoose(n, k)) +
                 pbeta(.5, s1, s2, log = TRUE) +
                 lbeta(s1, s2) - log(.5))

  # TTB is best model:
  expect_named(margs <- strategy_postprob(k, n, preds), names(preds))
  expect_gt(margs[1], max(margs[-c(1)]))

  # WADD is best model:
  expect_named(margs <- strategy_postprob(c(0,n[3],0), n, preds), names(preds))
  expect_gt(margs[2], max(margs[-c(2)]))

  # errors
  expect_error(strategy_marginal(k_false, n, multinomineq:::as_strategy(TTB)))
})

# correct results for first 5 participants
pp_heck2017 <-
  structure(c(0.0143891898881328, 0.00295077454126497, 0.0331602104316971,
              0.433027100317965, 0.000113503447049151, 0.55600835835707, 0.972993444432449,
              1.90597125805024e-07, 7.22294721347269e-15, 0.99606712429216,
              0.429602450108647, 0.0240556457350732, 0.966837178517799, 0.558558075453203,
              0.00381937225918344, 9.19234803419218e-26, 3.15727497160526e-30,
              2.38246407255183e-06, 0.00826643338052809, 5.58105937138287e-36,
              3.85184612683762e-28, 1.23933685130642e-32, 3.76238381306751e-08,
              0.000148390848296811, 2.2640821702952e-38, 1.6461504545141e-09,
              1.35291212484574e-07, 3.65467847146544e-10, 8.24139454860906e-22,
              1.60719296224333e-12, 1.42790147074372e-35, 4.39666823466277e-38,
              8.83421914025033e-18, 1.84856793333073e-27, 8.77819376577097e-45
  ), .Dim = c(5L, 7L), .Dimnames = list(c("101", "102", "103",
                                          "104", "105"), c("baseline", "WADD", "WADDprob", "TTB", "TTBprob",
                                                           "EQW", "GUESS")))

test_that("results match for Heck et al. (2017)", {
  data(heck2017)
  head(heck2017)
  n <- rep(40, 4)

  # cue validities and values
  v <- c(.9, .8, .7, .6)
  cueA <- matrix(c(-1,  1,  1, -1,
                   1, -1, -1,  1,
                   -1,  1,  1,  1,
                   1, -1, -1, -1),
                 ncol = 4, byrow = TRUE)
  cueB <- matrix(c(-1, -1, -1, -1,
                   -1, 1 , -1, 1 ,
                   -1, 1 , 1 , -1,
                   -1, 1 , 1 , -1),
                 ncol = 4, byrow = TRUE)

  # get predictions
  strategies <- c("baseline",  "WADD", "WADDprob",
                 "TTB",  "TTBprob", "EQW", "GUESS")
  strats <- strategy_multiattribute(cueA, cueB, v, strategies)
  pp <- unname(strategy_postprob(heck2017[1:5,], rep(40, 4), strats))
  heck <- unname(pp_heck2017)
  expect_equivalent(pp, heck)
})
