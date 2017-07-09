# vectors of predictions for 3 item types (Hilbig & Moshagen, 2014)
TTB <- c(-1, -1, -1)
WADD <- c(-1, 1, -1)
WADDprob <- c(-1, 3, -2)
EQW <- c(-1, 1, 0)
GUESS <- rep(0, 3)
baseline <- 1:3



k <- c(0, 1, 1)
k_false <- c(11,1,1)
n <- c(10, 15, 12)

test_that('ML estimation works', {

  expect_length(estimate_error(k, n, GUESS), 0)
  expect_equivalent(estimate_error(k, n, baseline), 1 - k/n)


  expect_silent(est_ttb <- maximize_ll(k, n, as_strategy(TTB)))
  expect_equal(est_ttb$error, sum(k)/sum(n))

  expect_silent(est_guess <- maximize_ll(k, n, as_strategy(GUESS)))
  expect_length(est_guess$error, 0)
  expect_equal(est_guess$loglik, sum(dbinom(k, n, .5, log = TRUE)))

  expect_silent(est_waddp <- maximize_ll(k, n, as_strategy(WADDprob)))
})

test_that('c_NML computation works', {
  n <- rep(2,3)
  # GUESS
  cnml_g <- compute_cnml(as_strategy(GUESS), n)
  expect_length(cnml_g, 6)
  expect_named(cnml_g, c("pattern","c", "ordered", "prior", "cnml", "n"))
  expect_equivalent(cnml_g$cnml, 0)

  # baseline + deterministic
  cnml_b <- compute_cnml(as_strategy(baseline, c=1, ordered = FALSE), n, c = 1)
  cnml_w <- compute_cnml(as_strategy(WADD), n)
  cnml_t <- compute_cnml(as_strategy(TTB), n)
  cnml_e <- compute_cnml(as_strategy(EQW), n)
  expect_equivalent(cnml_w$cnml, cnml_t$cnml)

  # probabilistic & NML order
  cnml_wp <- compute_cnml(as_strategy(WADDprob), n)
  expect_gt(cnml_b$cnml, cnml_wp$cnml)
  expect_gt(cnml_wp$cnml, cnml_w$cnml)
  expect_gt(cnml_w$cnml, cnml_e$cnml)
  expect_gt(cnml_e$cnml, cnml_g$cnml)
})

test_that('complexity for NML computation works', {
  n <- rep(2,3)
  # GUESS
  cnml_g <- compute_cnml(as_strategy(GUESS, prior = c(1.5, 1.5)), n)
  expect_length(cnml_g, 6)
  expect_equivalent(cnml_g$cnml, 0)

  # baseline + deterministic
  cnml_b <- compute_cnml(as_strategy(baseline, c = 1, ordered =FALSE, prior = c(1.5, 1.5)), n)
  cnml_w <- compute_cnml(as_strategy(WADD, prior = c(1.5, 1.5)), n)
  cnml_t <- compute_cnml(as_strategy(TTB, prior = c(1.5, 1.5)), n)
  cnml_e <- compute_cnml(as_strategy(EQW, prior = c(1.5, 1.5)), n)
  expect_equivalent(cnml_w$cnml, cnml_t$cnml)

  # probabilistic & NML order
  cnml_wp <- compute_cnml(as_strategy(WADDprob, prior = c(1.5, 1.5)), n)
  expect_gt(cnml_b$cnml, cnml_wp$cnml)
  expect_gt(cnml_wp$cnml, cnml_w$cnml)
  expect_gt(cnml_w$cnml, cnml_e$cnml)
  expect_gt(cnml_e$cnml, cnml_g$cnml)
})

