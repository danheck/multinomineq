#' Maximum-likelihood Estimate
#'
#' Get ML estimate for product-binomial/multinomial model with linear inequality constraints.
#'
#' @inheritParams count_binom
#' @inheritParams strategy_marginal
#' @param n.fit number of calls to \link[stats]{constrOptim}.
#' @param ... further arguments passed to the \code{control} list of \code{\link[stats]{constrOptim}} (e.g., \code{maxit = 5000}.
#' @return the list returned by the optimizer \code{\link[stats]{constrOptim}}.
#' @examples
#' # linear order: p1 < p2 < p3 < .50
#' # (cf. WADDprob in ?strategy_multiattribute)
#' A <- matrix(c(1, -1,  0,
#'               0,  1, -1,
#'               0,  0,  1),
#'             ncol = 3, byrow = TRUE)
#' b <- c(0, 0, .50)
#' ml_binom(k = c(4,1,23), n = 40, A, b)[1:2]
#' ml_multinom(k = c(4,36,  1,39,  23,17),
#'             options = c(2,2,2), A, b)[1:2]
#'
#'
#'
#' # probabilistic strategy:  A,A,A,B  [e1<e2<e3<e4<.50]
#' strat <- list(pattern = c(-1, -2, -3, 4),
#'               c = .5, ordered = TRUE, prior = c(1,1))
#' ml_binom(c(7,3,1, 19), 20, strategy = strat)[1:2]
#' @export
ml_binom <- function(k, n, A, b, map, strategy, n.fit = 1, start, ...){
  if (length(n) == 1) n <- rep(n, length(k))

  if (!missing(strategy)){
    pt <- strategy_to_Ab(strategy)
    A <- pt$A
    b <- pt$b
    map <- strategy$pattern
  } else {
    check_Abknprior(A, b, k, n)
    if (missing(map))
      map <- 1:ncol(A)
  }
  aggr <- map_k_to_A(k, n, A, map)
  options <- rep(2, ncol(A))
  k_all <- add_fixed(aggr$k, options = options, sum = aggr$n)
  ml_multinom(k_all, options, A, b, n.fit = n.fit, start = start,...)
}

#' @inheritParams count_multinom
#' @rdname ml_binom
#' @export
ml_multinom <- function(k, options, A, b, n.fit = 1, start, ...){
  check_Abokprior(A, b, options, k)
  # est <- k / c(tapply(k, rep(1:length(options), options), sum))
  # start <- drop_fixed(est, options)
  tmp <- Ab_multinom(options, A, b, nonneg = TRUE)
  A <- tmp$A
  b <- tmp$b
  if (missing(start) || !all(A %*% start < b))
    start <- find_inside(A, b, options = options)
  check_start(start, A, b, interior = TRUE)
  tryCatch (oo <- constrOptim(start, f = loglik_multinom, grad = grad_multinom,
                              k = k, options = options,
                              ui = - A, ci = - b, control = list(fnscale = -1, ...)),
            error = function(e) {print(e);
              cat("\n\n  (optimization failed with start values =", start, ")\n")})
  cnt <- 1
  while (cnt < n.fit){
    start <- find_inside(A, b, random = TRUE)
    oo2 <- constrOptim(start, loglik_multinom, grad = grad_multinom,
                       k = k, options = options, ui = - A, ci = - b,
                       control = list(fnscale = -1))
    cnt <- cnt + 1
    if (oo2$value > oo$value) oo <- oo2
  }
  names(oo$par) <- colnames(A)
  oo
}

# ' Choice Probabilities Implied by Strategy Predictions
# '
# ' Returns a vector with probabilities of choosing Option B for a predicted
# ' pattern and corresponding error probabilities.
# '
# ' @param error parameter vector of error probabilities (depending on \code{pattern}).
# '   A single value for deterministic strategies and an
# '   order vector of probabilities for probabilistic strategies (cf. examples).
# ' @inheritParams compute_marginal
# ' @examples
# ' # Predicted:  B,B,[guess],A
# ' # (deterministic: with constant error = .10)
# ' s1 <- list(pattern = c(1, 1, 0, -1),
# '            ordered = FALSE, c = .1, prior = c(1,1))
# ' error_to_prob(error = .10, s1)
# '
# ' # Predicted:  B,B,A,[guess]
# ' # (probabilistic: with ordered errors e1<e2<e3<.50)
# ' s2 <- list(pattern = c(1, 2, -3, 0),
# '            ordered = TRUE, c = .5, prior = c(1,1))
# ' error_to_prob(error = c(.10, .15,.20), s2)
# ' @export
error_to_prob <- function(error, strategy){
  pattern <- strategy$pattern
  if (all(pattern == 0))
    return(rep(.5, length(pattern)))
  if (!is.null(error) && any(error > strategy$c, error < 0))
    return(rep(NA, length(pattern)))
  # error per item type:
  if (strategy$ordered){
    if (any(sort(error) != error))
      return(rep(NA, length(pattern)))
    error_label <- rev(get_error_unique(pattern))
    error_per_type <- error[match(abs(pattern), error_label)]
  } else {
    error_per_type <- error
  }
  # reversed items and guessing:
  probB <- ifelse(pattern < 0, error_per_type,
                  ifelse(pattern > 0, 1 - error_per_type, .5))
  probB
}


# unique error indices/labels
get_error_unique <- function(pattern){
  pred <- pattern[pattern != 0]
  sort(unique(abs(pred)))
}

get_error_idx <- function (pattern){
  n_error <- get_error_number(pattern)
  e_idx <- rank(-abs(pattern), ties.method = "min")
  e_idx <- match(e_idx, sort(unique(e_idx[pattern != 0])))
  e_idx[pattern == 0] <- .50
  e_idx
}

get_error_number <- function(pattern){
  length(get_error_unique(pattern))
}

# product binomial loglikelihood with luckiness
#
# luckiness: -log(a(prob))   not normalized
# => complexity of order-constrained models can be compared
# default:    standard NML; luck=c(1,1)
# uniform BF: inverse of Jeffreys prior; luck=c(1.5, 1.5)
loglik <- function (error, k, n, strategy){
  pb <- error_to_prob(error, strategy)
  loglik_binom(pb, k, n)
}

loglik_binom <- function (p, k, n){
  try(ll <- sum(dbinom(x = k, size = n, prob = p, log = TRUE)), silent = TRUE)
  if (is.null(ll) || !is.na(ll) && ll == - Inf && all(k <= n))
    ll <- MIN_LL
  ll
}

loglik_multinom <- function (p, k, options){
  p_all <- add_fixed(p, options = options, sum = 1)
  ll <- sum(k * log(p_all))
  # tapply(p, oo, function(pp) dmultinom())
  # ll <- sum(dmultinom(x = k, size = n, prob = p, log = TRUE))
  if (!is.na(ll) && ll == - Inf)
    ll <- MIN_LL
  ll
}


grad_multinom <- function (p, k, options){
  p_all <- add_fixed(p, options = options, sum = 1)
  idx_K <- cumsum(options)
  K <- k[idx_K]
  pK <- p_all[idx_K]
  g <- k[-idx_K] / p - rep(K / pK, options - 1)
  # check: print(rbind(g = g, numderiv = numDeriv::grad(stratsel:::loglik_multinom, p, k=k, options=options)))
  # ll <- sum(dmultinom(x = k, size = n, prob = p, log = TRUE))
  # if (!is.na(g) && ll == - Inf)
  # g[!is.na(g) & g == ] <- MIN_LL
  g
}


# count adherence/error rates
# used to get ML estimate : adherence/n
#
# do not use strategy as input: priors for BF are counted elsewhere (clearer)
count_errors <- function (k, n, pattern, luck = c(1,1), ordered = FALSE, prob = TRUE){
  # reversed items:
  tmp <- k
  pred_B <- pattern > 0
  tmp[pred_B] <- n[pred_B] - k[pred_B]

  error_label <- get_error_unique(pattern)
  if (ordered)
    error_label <- rev(error_label)
  idx <- match(abs(pattern), error_label, nomatch = NA)
  cnt_predicted <- tapply(tmp, list(idx), sum) + luck[1] - 1
  if (prob){
    cnt_total <-  tapply(n, list(idx), sum) + sum(luck) - 2
    cnt_predicted <- cnt_predicted / cnt_total
  }
  unname(c(cnt_predicted), force = TRUE)
}


# move analytical estimate into interior of parameter space
adjust_error <- function (error, strategy, bound = 1e-10){
  error_adj <- error
  if (strategy$ordered){
    # simple linear order constraints
    error_adj <- adj_iterative(error, strategy$c, bound)
    error_adj <- sort(error_adj + runif(length(error), 0, bound/4))

    #### complex order constraints (ui, ci)
    # stop("complex order constraints (A, b) not yet implemented")
  } else if (length(error) > 0) {
    error_adj <- pmax(bound, pmin(strategy$c - bound, error))
  }
  error_adj
}

