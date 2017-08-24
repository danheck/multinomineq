#########################
# Compute NML complexity for decision strategies
#########################

# ' Maximum-likelihood Estimate
# '
# ' Get ML estimate (possible weighted by luckinesss function).
# '
# ' @inheritParams strategy_cnml
# ' @param k vector with observed frequencies for Option B
# ' @template details_strategy
# ' @examples
# ' k <- c(3,5,19, 9)  # frequency Option B
# ' n <- rep(20, 4)    # number of choices per type
# '
# ' # deterministic strategy: A,A,B,guess  [e<.50]
# ' s1 <- list(pattern = c(-1, -1, 1, 0),
# '            c = .5, ordered = FALSE, prior = c(1,1))
# ' maximize_ll(k, n, s1)  # = (3+5+20-19)/(3*20)
# '
# ' # probabilistic strategy:  A,A,A,B  [e1<e2<e3<e4<.50]
# ' s2 <- list(pattern = c(-1, -2, -3, 4),
# '            c = .5, ordered = TRUE, prior = c(1,1))
# ' maximize_ll(k, n, s2)
# ' @export
maximize_ll <- function(k, n, strategy, n.fit = 5){
  check_strategy(strategy)

  # analytic ML estimate, boundary correction and order constraints
  est <- count_errors(k, n, strategy$pattern, luck = strategy$prior,
                      ordered = strategy$ordered)
  start <- adjust_error(est, strategy)
  n_error <- length(est)
  # check for unordered models (EQW, TTB, WADD, baseline)
  # oo <- optim(ll_luckiness, lower=0, upper=c, control = list(fn=-1),
  #                k = k, n = n, pattern=pattern, luck = luck)

  if (n_error > 1 && strategy$ordered){
    start <- adjust_error(est, strategy, BOUND)
    # pt <- as_polytope(strategy)
    # ui <- - pt$A
    # ci <- - pt$b
    ui <- rbind(diag(n_error), 0) - rbind(0, diag(n_error))
    ci <- c(rep(0, n_error), - strategy$c)
    tryCatch (oo <- constrOptim(start, loglik_nml, grad = NULL,
                                k = k, n = n, strategy = strategy,
                                ui = ui, ci = ci, control = list(fnscale = -1)),
              error = function(e) {print(e);
                cat("\n\n  (optimization failed with start values =", start, ")\n")})
    cnt <- 1
    while (cnt < n.fit){
      start <- sort(runif(n_error, 0, strategy$c))
      oo2 <- constrOptim(start, loglik_nml, grad = NULL,
                         k = k, n = n, strategy = strategy,
                         ui = ui, ci = ci, control = list(fnscale = -1))
      oo2
      cnt <- cnt + 1
      if (oo2$value > oo$value) oo <- oo2
    }
    start <- oo$par
  }
  ll <- loglik_nml(start, k, n, strategy)
  list("loglik" = ll, "error" = start, "strategy" = strategy)
}


loglik_nml <- function (error, k, n, strategy){
  pb <- error_to_prob(error, strategy)
  ll <- sum(dbinom(x = k, size = n, prob = pb, log = TRUE))    # likelihood
  if (!is.null(error))
    ll <- ll + sum(dbeta(error, strategy$prior[1], strategy$prior[2], log = TRUE))  # luckiness
  if (!is.na(ll) && ll == - Inf && all(k <= n))
    ll <- MIN_LL
  ll
}

# Compute NML Complexity
#' @rdname strategy_nml
#' @export
strategy_cnml <- function(strategy, n, n.fit = 3, cores = 1){
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


#' Strategy Classification with Normalized Maximum Likelihood (NML)
#'
#' Computes the normalized maximum likelihood (NML), a model-selection method
#' based on the minimum description length principle. Essentially, NML penalizes
#' model complexity by summing the maximum likelihood values for all possible data vectors.
#'
# ' @inheritParams maximize_ll
#' @inheritParams strategy_marginal
#' @inheritParams strategy_postprob
#' @param n.fit number of repeated fitting runs with random starting values
#' @param cores number of processing units to be used
#' @param strategies list with strategies including the pre-computed
#'           NML complexity values (by using \code{strategy_cnml}).
#'           Note that the sample size \code{n} must match.
#' @details
#' Note that the prior parameters in \code{strategy} actually define the NML luckiness function.
#' (i.e., parameters of a beta distribution).
#' For luckiness NML, \code{prior = c(1.5, 1.5)} is asymptotically equivalent to
#' a uniform prior \code{prior = c(1, 1)} for the Bayes factor
#'
#' @template details_strategy
#' @examples
#' ### strategy predictions
#' # strategy: A/A/B with error probabilities {e2,e3,e1<.20}
#' s1 <- list(pattern = c(-3, -1, +2), c = .2,
#'            ordered = FALSE, prior = c(1,1))
#' # predict: A/A/B with ordered errors {e2<e3<e1<.50}
#' s2 <- list(pattern = c(-3, -1, +2), c = .5,
#'            ordered = TRUE, prior = c(1,1))
#'
#' ### compute NML complexity terms (depends on samples size!)
#' n <- c(2, 1, 2)
#' c1 <- strategy_cnml(s1, n)
#' c2 <- strategy_cnml(s2, n)
#'
#' ### NML model selection
#' nml <- strategy_nml(c(0,1,1), n, list(c1, c2))
#' model_weights(-nml)
#' @export
strategy_nml <- function(k, n, strategies, n.fit = 5, cores = 1){

  if (is.matrix(k) || is.data.frame(k)){
    if (cores > 1){
      cl <- makeCluster(cores)
      nml <- parApply(cl = cl, k, 1, strategy_nml, n = n,
                      strategies = strategies, n.fit = n.fit)
      stopCluster(cl)
    } else {
      nml <- apply(k, 1, strategy_nml, n = n,
                   strategies = strategies, n.fit = n.fit)
    }
    if (is.matrix(nml)) nml <- t(nml)

  } else {
    k <- unlist(k)
    n <- unlist(n)
    check_knp(k, n, n)
    sapply(strategies, check_cnml, n = n)
    lls <- sapply(strategies, function(ss)
      maximize_ll(k = k, n = n, ss, n.fit = n.fit)$loglik)

    nml <- - lls + sapply(strategies, "[[", "cnml")
    strat_names <- sapply(strategies, function(ss) attr(ss, "label"))
    sel <- !sapply(strat_names, is.null)
    names(nml)[sel] <- strat_names[sel]
  }
  nml
}


# ' ## (data)
# ' \dontrun{
# ' # strategy classification by NML (can take hours)
# ' cnmls <- strategy_cnml(preds, n, cores = 3)
# ' strategy_nml(heck2017[1:4,], n, cnmls)
# ' }
