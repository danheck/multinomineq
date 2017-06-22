
cueA
v <- seq(.9, .6, -.1)

TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
attr(WADDprob, "ordered") <- TRUE

b <- c(0, 1, 1)
b_false <- c(11,1,1)
n <- c(10, 15, 12)

test_that('predictions work as expected', {

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
})
