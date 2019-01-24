


test_that('multinomineq:::loglik_binom/_multinom gives correct results', {

  n <- rep(10,3)
  k <- c(3,2,5)
  expect_equal(multinomineq:::loglik_binom(c(.5,.5,.5), k ,n), - sum(n * log(.5)))
  p <- runif(3)
  expect_equal(multinomineq:::loglik_binom(p, k ,n), - sum(k * log(p), (n-k)*log(1-p)))
  expect_equal(suppressWarnings(multinomineq:::loglik_binom(c(0,0,0), k ,n)), Inf)

  p <- rpdirichlet(1, rep(1,7), c(3,4))
  p_all <- c(add_fixed(p, c(3,4)))
  k <- rpmultinom(p, 1, c(3,4))
  expect_equal(multinomineq:::loglik_multinom(p, c(k), c(3,4)), - sum(k * log(p_all)))
  expect_equal(suppressWarnings(multinomineq:::loglik_multinom(c(0,0,0,0,0), k, c(3,4))), Inf)
})


# vectors of predictions for 3 item types (Hilbig & Moshagen, 2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-3, 1, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3

k <- c(0, 5, 1)
k_false <- c(11,1,1)
n <- c(10, 15, 12)

test_that('predictions work as expected', {

  # deterministic models
  expect_equal(multinomineq:::get_error_unique(TTB), 1)
  expect_length(multinomineq:::get_error_unique(GUESS), 0)
  # expect_equal(multinomineq:::error_number(TTB), 1)
  expect_equal(multinomineq::: error_to_prob(.123, multinomineq:::as_strategy(TTB)), rep(.123, 3))

  # probabilistic models
  expect_equal(multinomineq:::get_error_unique(WADDprob), 1:3)
  # expect_equal(multinomineq:::error_number(WADDprob), 3)
  expect_equal(multinomineq::: error_to_prob(c(.1,.25, .3), multinomineq:::as_strategy(WADDprob)),
               c(.1, 1-.3, .25))
  # baseline
  expect_equal(multinomineq::: error_to_prob(c(.1,.25, .3), multinomineq:::as_strategy(baseline, c=1, ordered=FALSE)),
               1-c(.1, .25, .3))
})

test_that("multinomineq::loglik_strategy works", {

  # loglik_strategy is now DEPRECATED [ was used for NML]
  expect_gte(multinomineq:::loglik_strategy(runif(1,0.1,.5), k ,n, multinomineq:::as_strategy(TTB)), 0)
  expect_lte(multinomineq:::loglik_strategy(sum(k)/sum(n), k ,n, multinomineq:::as_strategy(TTB)),
             multinomineq:::loglik_strategy(runif(1,0,.5), k ,n, multinomineq:::as_strategy(TTB)))

  expect_gte(multinomineq:::loglik_strategy(runif(1,0,.5), k ,n, multinomineq:::as_strategy(WADD)), 0)
  expect_gte(multinomineq:::loglik_strategy(runif(3,0,1), k ,n,
                                   multinomineq:::as_strategy(baseline, c=1, ordered = FALSE)), 0)
  expect_gte(multinomineq:::loglik_strategy(c(), k ,n, multinomineq:::as_strategy(GUESS)), 0)
  expect_equal(multinomineq:::loglik_strategy(c(), k ,n, multinomineq:::as_strategy(GUESS)),
               - sum(n * log(.5)))

  # expected errors
  expect_identical(suppressWarnings(
    multinomineq:::loglik_strategy(0, k_false , n, multinomineq:::as_strategy(TTB))), Inf)
  # inadmissible parameters
  expect_identical(suppressWarnings(
    multinomineq:::loglik_strategy(1.5, k , n, multinomineq:::as_strategy(TTB))), NA_real_)
  expect_identical(suppressWarnings(
    multinomineq:::loglik_strategy(-.1, k , n, multinomineq:::as_strategy(TTB))), NA_real_)
  expect_identical(suppressWarnings(
    multinomineq:::loglik_strategy(-.1, k , n, multinomineq:::as_strategy(WADDprob))), NA_real_)
  expect_identical(suppressWarnings(
    multinomineq:::loglik_strategy(c(.5,.4,.3), k , n, multinomineq:::as_strategy(WADDprob))), NA_real_)

})

