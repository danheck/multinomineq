library("testthat")
set.seed(123)

# models of Hilbig & Moshagen (2014)
v <- c(.9, .8, .7, .6)
cuesA <- matrix(c(1, 1, 1, -1,
                  1, -1, -1, -1,
                  1, 1, 1, -1), ncol = 4, byrow = TRUE)
cuesB <- matrix(c(-1,1,-1,1,
                  -1,1,1,-1,
                  -1,1,1,1), ncol = 4, byrow = TRUE)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-2.63905733,  0.03636764 ,-1.79175947)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
predictions <- list("baseline" = baseline, "WADDprob" = WADDprob,
                    "TTB" = TTB, "WADD" = WADD,
                    "EQW" = EQW, "GUESS" = GUESS)
strats <- names(predictions)

set.seed(12345)

test_that('predictions work for one item type', {

  ca <- c(1, 1, -1, 1)
  cb <- c(-1, -1, -1, -1)
  expect_named(p1 <- strategy_multiattribute(ca, cb, v, "TTBprob"),
               c("pattern", "c"   ,    "ordered", "prior",   "label"))
  expect_equivalent(p1$pattern, -1)
  expect_equivalent(strategy_multiattribute(ca, cb, v, "WADDprob")$pattern,
                    sum(log( 1/c(.9,.8,.6)-1 )))
  expect_equal(strategy_multiattribute(ca, cb, v, "WADD")$pattern, -1)
  expect_equal(strategy_multiattribute(ca, cb, v, "EQW")$pattern, -1)
  expect_equal(strategy_multiattribute(ca, cb, v, "GUESS")$pattern, 0)
  expect_named(strategy_multiattribute(ca, cb, v, strats), strats)

  # deterministic models
  expect_equal(multinomineq:::get_error_unique(TTB), 1)
  expect_length(multinomineq:::get_error_unique(GUESS), 0)
  expect_equal(multinomineq:::get_error_number(TTB), 1)
  expect_equal(multinomineq:::error_to_prob(.123, multinomineq:::as_strategy(TTB)), rep(.123, 3))

  # probabilistic models
  expect_equal(multinomineq:::get_error_unique(WADDprob), sort(abs(WADDprob)))
  expect_equal(multinomineq:::get_error_number(WADDprob), 3)
  expect_equal(multinomineq:::error_to_prob(c(.1,.2,.45), multinomineq:::as_strategy(WADDprob)),
               c(.1, 1 - .45, .2))
  # baseline
  e <- c(.1 , .25, .3)
  expect_equal(multinomineq:::error_to_prob(e, multinomineq:::as_strategy(baseline, ordered = FALSE, c=1)), 1 - e)

})

test_that('predictions work for multiple item types', {
  expect_named(p <- strategy_multiattribute(cuesA, cuesB, v, strats), strats)
  expect_equal(lapply(p, "[[", "pattern"), predictions)

  data(heck2017_raw)
  cA <- heck2017_raw[,paste0("a",1:4)]
  cB <- heck2017_raw[,paste0("b",1:4)]
  v <- c(.9, .8, .7, .6)
  ttb <- strategy_multiattribute(cA, cB, v, "TTB")
  expect_equivalent(heck2017_raw$ttb , ifelse(ttb$pattern == -1, "A", "B"))
  wadd <- strategy_multiattribute(cA, cB, v, "WADD")
  expect_equivalent(heck2017_raw$wadd , ifelse(wadd$pattern == -1, "A", "B"))
  waddp <- strategy_multiattribute(cA, cB, v, "WADDprob")
  expect_equivalent(heck2017_raw$logoddsdiff , - waddp$pattern)
  eqw <- strategy_multiattribute(cA, cB, v, "EQW")
  expect_equivalent(heck2017_raw$eqw , ifelse(eqw$pattern == -1, "A",
                                              ifelse(eqw$pattern == 1, "B", "guess")))
  ttbp <- strategy_multiattribute(cA, cB, v, "TTBprob")
  expect_equivalent(unname(1/ heck2017_raw$ttbsteps),
                    unname(abs(ttbp$pattern)))
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

test_that("predictions/frequency table matches for Heck et al. (2017)", {

  cueA <- heck2017_raw[,paste0("a",1:4)]
  cueB <- heck2017_raw[,paste0("b",1:4)]
  v <- c(.9, .8, .7, .6)

  strat_labels <- c("TTB", "TTBprob", "WADD",
                    "WADDprob", "EQW", "GUESS")
  strat <- strategy_multiattribute(cueA, cueB, v, strat_labels)
  types <- strategy_unique(strat)

  item_rev <- paste(heck2017_raw$itemtype,
                    heck2017_raw$reversedorder)
  expect_equal(unique(table(item_rev, types$item_type)), c(0, 2080))

  freq <- with(heck2017_raw,
               table(vp, types$item_type, choice))
  freqB <- freq[,4:1,1] + # reversed items: Option A
    freq[,5:8,2]   # non-rev. items: Option B
  expect_true( all( (40 - freqB)[,c(3,2,4,1)] - heck2017 == 0))

  pp <- strategy_postprob(freqB[1:5,], rep(40, 4),
                          types$strategies)
  expect_equivalent(pp[,strat_labels], pp_heck2017[,strat_labels])
})


test_that("strategy_to_Ab returns the correct A/b representation",{
  # compare results to encompassing BF method:
  encompassing <- list(pattern = 1:4, c = 1, ordered = FALSE, prior = c(1,1))
  k <- c(5, 4, 2, 0)
  n <- rep(10, 4)

  strat <- list(pattern = -c(1:4),  # A,A,A,A  e4<e3<e2<e1<.5  | .5>t1>t2>t3>t4
                c = .5, ordered = TRUE,
                prior = c(1,1))
  pt <- strategy_to_Ab(strat)
  m1 <- strategy_postprob(k, n, list(strat, encompassing))

  bf <- bf_binom(k, n, pt$A, pt$b, M = 1e6, log = TRUE)
  expect_equal(log(m1[1] / m1[2]), bf["bf_0u",1], tol = 3*bf["bf_0u",2])

  strat <- list(pattern =-c(4:1),  # A,A,A,A  e1<e2<e3<e4<.5
                c = .5, ordered = TRUE,
                prior = c(1,1))
  pt <- strategy_to_Ab(strat)
  m1 <- strategy_postprob(k, n, list(strat, encompassing))
  bf <- bf_binom(k, n, pt$A, pt$b, M = 1e6, log = TRUE)
  expect_equal(log(m1[1] / m1[2]), bf["bf_0u",1], tol = 3*bf["bf_0u",2])

  strat <- list(pattern =c(1,-5,2,-3),  # A,A,A,A  e1<e2<e3<e4<.5
                c = .5, ordered = TRUE,
                prior = c(1,1))
  pt <- strategy_to_Ab(strat)
  m1 <- strategy_postprob(k, n, list(strat, encompassing))
  bf <- bf_binom(k, n, pt$A, pt$b, M = 1e6, log = TRUE)
  expect_equal(log(m1[1] / m1[2]), bf["bf_0u",1], tol = 3*bf["bf_0u",2])
})
