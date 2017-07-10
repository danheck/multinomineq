# vectors of predictions for 3 item types (Hilbig & Moshagen, 2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-3, 1, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3

k <- c(0, 1, 1)
k_false <- c(11,1,1)
n <- c(10, 15, 12)

test_that('predictions work as expected', {

  # deterministic models
  expect_equal(get_error_unique(TTB), 1)
  expect_length(get_error_unique(GUESS), 0)
  expect_equal(get_error_number(TTB), 1)
  expect_equal(error_to_prob(.123, as_strategy(TTB)), rep(.123, 3))

  # probabilistic models
  expect_equal(get_error_unique(WADDprob), 1:3)
  expect_equal(get_error_number(WADDprob), 3)
  expect_equal(error_to_prob(c(.1,.25, .3), as_strategy(WADDprob)),
               c(.1, 1-.3, .25))
  # baseline
  expect_equal(error_to_prob(c(.1,.25, .3), as_strategy(baseline, c=1, ordered=FALSE)),
               1-c(.1, .25, .3))
})

test_that('loglik works', {

  expect_lte(loglik(runif(1,0.1,.5), k ,n, as_strategy(TTB)), 0)
  expect_gte(loglik(sum(k)/sum(n), k ,n, as_strategy(TTB)),
             loglik(runif(1,0,.5), k ,n, as_strategy(TTB)))

  expect_lte(loglik(runif(1,0,.5), k ,n, as_strategy(WADD)), 0)
  expect_lte(loglik(runif(3,0,1), k ,n, as_strategy(baseline, c=1, ordered = FALSE)), 0)
  expect_lte(loglik(c(), k ,n, as_strategy(GUESS)), 0)

  # expected errors
  expect_identical(loglik(runif(1,0,.5), k_false , n, as_strategy(TTB)), -Inf)
  # inadmissible parameters
  expect_identical(loglik(1.5, k , n, as_strategy(TTB)), NA_real_)
  expect_identical(loglik(-.1, k , n, as_strategy(TTB)), NA_real_)
  expect_identical(loglik(-.1, k , n, as_strategy(WADDprob)), NA_real_)
  expect_identical(loglik(c(.5,.4,.3), k , n, as_strategy(WADDprob)), NA_real_)

  # standard ML estimation
  expect_silent(ml <- maximize_ll(k, n, as_strategy(rep(-1, length(k)))))
  est <- sum(k)/sum(n)
  expect_equal(ml$error, est)
  expect_equal(ml$loglik, sum(dbinom(k, n, est, log = TRUE)))

  # luckiness ML estimation
  expect_silent(ml <- maximize_ll(k, n, as_strategy(rep(-1, length(k)),
                                                    prior = c(1.5, 1.5)))) # Beta(1.5, 1.5)
  est <- (sum(k) + .5)/(sum(n) + 1)
  expect_equal(ml$error, est)
  expect_equal(ml$loglik, dbeta(est,1.5, 1.5, log = TRUE) +
                 sum(dbinom(k, n, est, log = TRUE)))

  # baseline luckiness
  expect_silent(ml <- maximize_ll(k, n, list(pattern = - seq_along(k), c = 1,
                                             ordered = FALSE, prior = c(1.5, 1.5))))
  est <- (k + .5)/(n + 1)
  expect_equal(ml$error, est)
  expect_equal(ml$loglik, sum(dbeta(est,1.5, 1.5, log = TRUE) +
                                dbinom(k, n, est, log = TRUE)))
})

