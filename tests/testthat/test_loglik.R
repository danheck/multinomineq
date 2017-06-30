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
  expect_equal(ml$loglik, sum(dbinom(k, n, est, log = TRUE)))

})

test_that('ML estimation works', {

  expect_length(estimate_par(k, n, GUESS), 0)
  expect_equivalent(estimate_par(k, n, baseline), 1 - k/n)


  expect_silent(est_ttb <- maximize_ll(k, n, TTB))
  expect_equal(est_ttb$est, sum(k)/sum(n))

  expect_silent(est_guess <- maximize_ll(k, n, GUESS))
  expect_length(est_guess$est, 0)
  expect_equal(est_guess$loglik, sum(dbinom(k, n, .5, log = TRUE)))

  expect_silent(est_waddp <- maximize_ll(k, n, WADDprob))


})

test_that('c_NML computation works', {
  n <- rep(2,3)
  # GUESS
  cnml_g <- compute_cnml(GUESS, n)
  expect_length(cnml_g, 6)
  expect_named(cnml_g, c("cnml", "n", "prediction", "c", "luck", "time"))
  expect_equal(cnml_g$cnml, 0)

  # baseline + deterministic
  cnml_b <- compute_cnml(baseline, n, c = 1)
  cnml_w <- compute_cnml(WADD, n)
  cnml_t <- compute_cnml(TTB, n)
  cnml_e <- compute_cnml(EQW, n)
  expect_equal(cnml_w$cnml, cnml_t$cnml)

  # probabilistic & NML order
  cnml_wp <- compute_cnml(WADDprob, n)
  expect_gt(cnml_b$cnml, cnml_wp$cnml)
  expect_gt(cnml_wp$cnml, cnml_w$cnml)
  expect_gt(cnml_w$cnml, cnml_e$cnml)
  expect_gt(cnml_e$cnml, cnml_g$cnml)

  # expect_equal(cnml_wp$cnml, 1.52748)
})

test_that('c_LNML computation works', {
  n <- rep(2,3)
  # GUESS
  cnml_g <- compute_cnml(GUESS, n, luck = c(1.5, 1.5))
  expect_length(cnml_g, 6)
  expect_equal(cnml_g$cnml, 0)

  # baseline + deterministic
  cnml_b <- compute_cnml(baseline, n, c = 1, luck = c(1.5, 1.5))
  cnml_w <- compute_cnml(WADD, n, luck = c(1.5, 1.5))
  cnml_t <- compute_cnml(TTB, n, luck = c(1.5, 1.5))
  cnml_e <- compute_cnml(EQW, n, luck = c(1.5, 1.5))
  expect_equal(cnml_w$cnml, cnml_t$cnml)

  # probabilistic & NML order
  cnml_wp <- compute_cnml(WADDprob, n, luck = c(1.5, 1.5))
  expect_gt(cnml_b$cnml, cnml_wp$cnml)
  expect_gt(cnml_wp$cnml, cnml_w$cnml)
  expect_gt(cnml_w$cnml, cnml_e$cnml)
  expect_gt(cnml_e$cnml, cnml_g$cnml)

  # expect_equal(cnml_wp$cnml, 1.52748)
})
