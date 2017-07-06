#########################
# Compute NML complexity for decision strategies
#########################

# -log(a(theta)); luckiness not normalized
# => complexity of order-constrained models can be compared
# default:    standard NML; luck=c(1,1)
# uniform BF: inverse of Jeffreys prior; luck=c(1.5, 1.5)
luckiness <- function(par, luck = c(1, 1)){
  dd <- 0
  if (length(par) > 0)
    dd <- sum(dbeta(par, luck[1], luck[2], log = TRUE))
  dd
}

ll_luckiness <- function(error, k, n, pattern, luck = c(1,1)){
  loglik(error, k, n, pattern) + luckiness(error, luck)
}

#' Maximum-likelihood Estimate
#'
#' Get ML estimate (possible weighted by luckinesss function).
#'
#' @inheritParams compute_cnml
#' @param k vector with observed frequencies for Option B
#' @template details_strategy
#' @examples
#' k <- c(3,5,19, 9)  # frequency Option B
#' n <- rep(20, 4)    # number of choices per type
#'
#' # deterministic strategy: A,A,B,guess  [e<.50]
#' s1 <- list(pattern = c(-1, -1, 1, 0),
#'            c = .5, ordered = FALSE, prior = c(1,1))
#' maximize_ll(k, n, s1)  # = (3+5+20-19)/(3*20)
#'
#' # probabilistic strategy:  A,A,A,B  [e1<e2<e3<e4<.50]
#' s2 <- list(pattern = c(-1, -2, -3, 4),
#'            c = .5, ordered = TRUE, prior = c(1,1))
#' maximize_ll(k, n, s2)
#' @export
maximize_ll <- function(k, n, strategy, n.fit = 5){
  check_strategy(strategy)
  luck <- strategy$prior
  pattern <- strategy$pattern
  c <- strategy$c

  # analytic ML estimate, boundary correction and order constraints
  est <- estimate_error(k, n, pattern, luck = luck)
  start <- adjust_error(est, pattern, c)
  n_error <- length(est)

  # check for unordered models (EQW, TTB, WADD, baseline)
  # oo <- optim(ll_luckiness, lower=0, upper=c, control = list(fn=-1),
  #                k = k, n = n, pattern=pattern, luck = luck)
  if (n_error > 1 && !is.null(attr(pattern, "ordered"))){
    start <- adjust_error(est, pattern, c, BOUND)
    # pt <- as_polytope(pattern, c)
    # ui <- - pt$A
    # ci <- - pt$b
    ui <- rbind(diag(n_error), 0) - rbind(0, diag(n_error))
    ci <- c(rep(0, n_error), -c)
    tryCatch (oo <- constrOptim(start, ll_luckiness, grad = NULL,
                                k = k, n = n, pattern = pattern, luck = luck,
                                ui = ui, ci = ci, control = list(fnscale = -1)),
              error = function(e) {print(e);
                cat("\n\n  (optimization failed with start values =", start, ")\n")})
    cnt <- 1
    while (cnt < n.fit){
      start <- sort(runif(n_error, 0, c))
      oo2 <- constrOptim(start, ll_luckiness, grad = NULL,
                         k = k, n = n, pattern = pattern, luck = luck,
                         ui = ui, ci = ci, control = list(fnscale = -1))
      oo2
      cnt <- cnt + 1
      if (oo2$value > oo$value) oo <- oo2
    }
    start <- oo$par
  }
  ll <- ll_luckiness(start, k, n, pattern, luck)
  list(loglik = ll, est = start, luck = luck)
}

#' Compute NML Complexity
#'
#' Enumerates discrete data space and sums the maximum likelihood values of all data vectors.
#'
#' @inheritParams compute_marginal
#' @param n.fit number of repeated fitting runs with random starting values
#' @param cores number of processing units to be used
#'
#' @details
#' Note that the prior parameters in \code{strategy} define the NML luckiness function.
#' (i.e., parameters of a beta distribution).
#' \code{luck = c(1.5, 1.5)} is equivalent to uniform prior for the Bayes factor
#'
#' @template details_strategy
#' @examples
#' # strategy: A/A/B with error probabilities {e2,e3,e1<.20}
#' s1 <- list(pattern = c(-3, -1, +2), c = .2,
#'            ordered = FALSE, prior = c(1,1))
#' compute_cnml(strat, n = c(2, 2, 2))
#'
#' # predict: A/A/B with ordered errors {e2<e3<e1<.50}
#' s2 <- list(pattern = c(-3, -1, +2), c = .5,
#'            ordered = TRUE, prior = c(1,1))
#' compute_cnml(s2, n = c(2, 2, 2))
#' @export
compute_cnml <- function(strategy, n, n.fit = 3, cores = 1){
  check_data_strategy(n, n, strategy)

  # generate all possible data and fit models
  dat <- do.call("expand.grid", lapply(n, function(nn) 0:nn))  # TODO: improve
  time <- system.time({
    if (cores > 1){
      cl <- makeCluster(cores)
      clusterExport(cl, c( "luck", "n.fit", "c", "n"),
                    envir=environment())
      lls <- parSapplyLB(cl, as.list(data.frame(t(dat))),
                         function(xx) maximize_ll(k = xx, n = n, strategy = strategy,
                                                  n.fit = n.fit)$loglik)
      stopCluster(cl)
    } else {
      lls <- apply(X = dat, MARGIN = 1,
                   function(xx) maximize_ll(k = xx, n = n, strategy = strategy,
                                            n.fit = n.fit)$loglik)
    }
  })

  strategy$cnml <- log(sum(exp(lls)))
  attr(strategy$cnml, "time") <- time["elapsed"]
  strategy$n <- n
  strategy
}


#' Strategy Selection Using NML
#'
#' Compute (L)NML for observed frequencies.
#'
#' @inheritParams compute_marginal
#' @inheritParams select_bf
#' @param strategy.list list with strategies including the pre-computed
#'           NML complexity values (see \code{\link{compute_cnml}}).
#'           Note that the sample size \code{n} must match.
#' @examples
#' n <- c(2, 2, 2)
#' # strategy: A/A/B with error probabilities {e2,e3,e1<.20}
#' s1 <- list(pattern = c(-3, -1, +2), c = .2,
#'            ordered = FALSE, prior = c(1,1))
#' s1n <- compute_cnml(s1, n = n)
#'
#' # predict: A/A/B with ordered errors {e2<e3<e1<.50}
#' s2 <- list(pattern = c(-3, -1, +2), c = .5,
#'            ordered = TRUE, prior = c(1,1))
#' s2n <- compute_cnml(s2, n = n)
#'
#' select_nml(c(0,2,2), n, list(s1n, s2n))
#' @export
select_nml <- function(k, n, strategy.list, n.fit = 5, cores = 1){

  if (is.matrix(k) || is.data.frame(k)){
    if (cores > 1){
      cl <- makeCluster(cores)
      nml <- parApply(cl = cl, k, 1, select_nml, n = n,
                      strategy.list = strategy.list, n.fit = n.fit)
      stopCluster(cl)
    } else {
      nml <- apply(k, 1, select_nml, n = n,
                   strategy.list = strategy.list, n.fit = n.fit)
    }
    if (is.matrix(nml)) nml <- t(nml)

  } else {
    k <- unlist(k)
    n <- unlist(n)
    check_knp(k, n, n)
    sapply(strategy.list, check_cnml, n = n)
    lls <- sapply(strategy.list, function(ss)
      maximize_ll(k = k, n = n, ss, n.fit = n.fit)$loglik)

    nml <- - lls + sapply(strategy.list, "[[", "cnml")
    strat_names <- sapply(strategy.list, function(ss) attr(ss, "label"))
    sel <- !sapply(strat_names, is.null)
    names(nml)[sel] <- strat_names[sel]
  }
  nml
}
