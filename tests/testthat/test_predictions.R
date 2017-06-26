# models of Hilbig & Moshagen (2014)
cuesA <- matrix(c(1, 1, 1, -1,
                  1, -1, -1, -1,
                  1, 1, 1, -1), ncol = 4, byrow = TRUE)
cuesB <- matrix(c(-1,1,-1,1,
                  -1,1,1,-1,
                  -1,1,1,1), ncol = 4, byrow = TRUE)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
attr(WADDprob, "ordered") <- TRUE


test_that('predictions work as expected', {

  # one item type
  ca <- c(1, 1, -1, 1)
  cb <- c(-1, -1, -1, -1)
  v <- c(.9, .8, .7, .6)
  expect_equivalent(predict_multiattribute(ca, cb, v, "TTBprob"), -1)
  expect_equivalent(predict_multiattribute(ca, cb, v, "WADDprob"),
                    sum(log( 1/c(.9,.8,.6)-1 )))
  expect_equal(predict_multiattribute(ca, cb, v, "WADD"), -1)
  expect_equal(predict_multiattribute(ca, cb, v, "EQW"), -1)
  expect_equal(predict_multiattribute(ca, cb, v, "GUESS"), 0)


  # deterministic models
  expect_equal(get_par_unique(TTB), 1)
  expect_length(get_par_unique(GUESS), 0)
  expect_equal(get_par_number(TTB), 1)
  expect_equal(get_prob_B(.123, TTB), rep(.123, 3))

  # probabilistic models
  expect_equal(get_par_unique(WADDprob), 1:3)
  expect_equal(get_par_number(WADDprob), 3)
  expect_equal(get_prob_B(c(.1,.25, .3), WADDprob),
               c(.1, .7, .25))
  # baseline
  expect_equal(get_prob_B(c(.1,.25, .3), baseline),
               1-c(.1, .25, .3))

  expect_equal(predict_multiattribute(cuesA, cuesB, v, "TTB"), TTB)
  expect_equal(predict_multiattribute(cuesA, cuesB, v, "WADD"), WADD)
  expect_equal(predict_multiattribute(cuesA, cuesB, v, "EQW"), EQW)
  expect_equal(predict_multiattribute(cuesA, cuesB, v, "GUESS"), GUESS)

  data(heck2017_raw)
  cA <- heck2017_raw[,paste0("a",1:4)]
  cB <- heck2017_raw[,paste0("b",1:4)]
  v <- c(.9, .8, .7, .6)
  ttb <- predict_multiattribute(cA, cB, v, "TTB")
  expect_equivalent(2*prop.table(table(heck2017_raw$ttb, ttb)),
                    as.table(diag(2)))
  wadd <- predict_multiattribute(cA, cB, v, "WADD")
  expect_equivalent(2*prop.table(table(heck2017_raw$wadd, wadd)),
                    as.table(diag(2)))
  waddp <- predict_multiattribute(cA, cB, v, "WADDprob")
  expect_equivalent(8*prop.table(table(heck2017_raw$logoddsdiff, -waddp)),
                    as.table(diag(8)))
  eqw <- predict_multiattribute(cA, cB, v, "EQW")
  expect_equal(sum(table(heck2017_raw$eqw, eqw) == 0), 6)

  ttbp <- predict_multiattribute(cA, cB, v, "TTBprob")
  expect_equal(sum(table(heck2017_raw$ttbsteps, ttbp) == 0), 6*2)

})
