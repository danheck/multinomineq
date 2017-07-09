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

test_that('predictions work for one item type', {

  ca <- c(1, 1, -1, 1)
  cb <- c(-1, -1, -1, -1)
  expect_named(p1 <- predict_multiattribute(ca, cb, v, "TTBprob"),
                c("pattern", "c"   ,    "ordered", "prior",   "label"))
  expect_equivalent(p1$pattern, -1)
  expect_equivalent(predict_multiattribute(ca, cb, v, "WADDprob")$pattern,
                    sum(log( 1/c(.9,.8,.6)-1 )))
  expect_equal(predict_multiattribute(ca, cb, v, "WADD")$pattern, -1)
  expect_equal(predict_multiattribute(ca, cb, v, "EQW")$pattern, -1)
  expect_equal(predict_multiattribute(ca, cb, v, "GUESS")$pattern, 0)
  expect_named(predict_multiattribute(ca, cb, v, strats), strats)

  # deterministic models
  expect_equal(get_error_unique(TTB), 1)
  expect_length(get_error_unique(GUESS), 0)
  expect_equal(get_error_number(TTB), 1)
  expect_equal(error_to_prob(.123, TTB), rep(.123, 3))

  # probabilistic models
  expect_equal(get_error_unique(WADDprob), sort(abs(WADDprob)))
  expect_equal(get_error_number(WADDprob), 3)
  expect_equal(error_to_prob(c(.1,.2,.45), as_strategy(WADDprob)),
               c(.1, 1 - .45, .2))
  # baseline
  e <- c(.1 , .25, .3)
  expect_equal(error_to_prob(e, as_strategy(baseline, ordered = FALSE, c=1)), 1 - e)

})

test_that('predictions work for multiple item types', {
  expect_named(p <- predict_multiattribute(cuesA, cuesB, v, strats), strats)
  expect_equal(lapply(p, "[[", "pattern"), predictions)

  data(heck2017_raw)
  cA <- heck2017_raw[,paste0("a",1:4)]
  cB <- heck2017_raw[,paste0("b",1:4)]
  v <- c(.9, .8, .7, .6)
  ttb <- predict_multiattribute(cA, cB, v, "TTB")
  expect_equivalent(heck2017_raw$ttb , ifelse(ttb$pattern == -1, "A", "B"))
  wadd <- predict_multiattribute(cA, cB, v, "WADD")
  expect_equivalent(heck2017_raw$wadd , ifelse(wadd$pattern == -1, "A", "B"))
  waddp <- predict_multiattribute(cA, cB, v, "WADDprob")
  expect_equivalent(heck2017_raw$logoddsdiff , - waddp$pattern)
  eqw <- predict_multiattribute(cA, cB, v, "EQW")
  expect_equivalent(heck2017_raw$eqw , ifelse(eqw$pattern == -1, "A",
                                              ifelse(eqw$pattern == 1, "B", "guess")))
  ttbp <- predict_multiattribute(cA, cB, v, "TTBprob")
  expect_equivalent(heck2017_raw$ttbsteps, abs(ttbp$pattern))

})
