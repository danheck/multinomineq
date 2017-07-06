# vectors of predictions for 3 item types (Hilbig & Moshagen, 2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3
attr(WADDprob, "ordered") <- TRUE

k <- c(0, 1, 1)
k_false <- c(11,1,1)
n <- c(10, 15, 12)

test_that('predictions work as expected', {

  # deterministic models
  expect_equal(get_error_unique(TTB), 1)
  expect_length(get_error_unique(GUESS), 0)
  expect_equal(get_error_number(TTB), 1)
  expect_equal(error_to_probB(.123, TTB), rep(.123, 3))

  # probabilistic models
  expect_equal(get_error_unique(WADDprob), 1:3)
  expect_equal(get_error_number(WADDprob), 3)
  expect_equal(error_to_probB(c(.1,.25, .3), WADDprob),
               c(.1, .7, .25))
  # baseline
  expect_equal(error_to_probB(c(.1,.25, .3), baseline),
               1-c(.1, .25, .3))
})

test_that('loglik works', {

  expect_lte(loglik(runif(1), k ,n, TTB), 0)
  expect_gte(loglik(sum(k)/sum(n), k ,n, TTB),
             loglik(runif(1), k ,n, TTB))

  expect_lte(loglik(runif(3), k ,n, WADD), 0)
  expect_lte(loglik(runif(3), k ,n, baseline), 0)
  expect_lte(loglik(c(), k ,n, GUESS), 0)

  # expected errors
  expect_identical(loglik(runif(1), k_false , n, TTB), -Inf)

  # standard ML estimation
  expect_silent(ml <- maximize_ll(k, n, rep(-1, length(k))))
  est <- sum(k)/sum(n)
  expect_equal(ml$est, est)
  expect_equal(ml$loglik, sum(dbinom(k, n, est, log = TRUE)))

  # luckiness ML estimation
  expect_silent(ml <- maximize_ll(k, n, rep(-1, length(k)),
                                  luck = c(1.5, 1.5)))      # Beta(1.5, 1.5)
  est <- (sum(k) + .5)/(sum(n) + 1)
  expect_equal(ml$est, est)
  expect_equal(ml$loglik, dbeta(est,1.5, 1.5, log = TRUE) +
                 sum(dbinom(k, n, est, log = TRUE)))

  # baseline luckiness
  expect_silent(ml <- maximize_ll(k, n, - seq_along(k),
                                  luck = c(1.5, 1.5), c = 1))
  est <- (k + .5)/(n + 1)
  expect_equal(ml$est, est)
  expect_equal(ml$loglik, sum(dbeta(est,1.5, 1.5, log = TRUE) +
                 dbinom(k, n, est, log = TRUE)))
})

