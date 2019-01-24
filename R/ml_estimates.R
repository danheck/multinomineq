#' Maximum-likelihood Estimate
#'
#' Get ML estimate for product-binomial/multinomial model with linear inequality constraints.
#'
#' @inheritParams count_binom
#' @inheritParams strategy_marginal
#' @param n.fit number of calls to \link[stats]{constrOptim}.
#' @param ... further arguments passed to the function
#'     \code{\link[stats]{constrOptim}}. To ensure high accuracy, the number of
#'     maximum iterations should be sufficiently large (e.g., by setting
#'     \code{control = list(maxit = 1e6, reltol=.Machine$double.eps^.6), outer.iterations = 1000}.
#'
#' @details
#' First, it is checked whether the unconstrained maximum-likelihood estimator
#' (e.g., for the binomial: \code{k/n}) is inside the constrained parameter space.
#' Only if this is not the case, nonlinear optimization with convex linear-inequality
#' constrained is used to estimate (A) the probability parameters \eqn{\theta}
#' for the Ab-representation or (B) the mixture weights \eqn{\alpha} for the V-representation.
#'
#' @return the list returned by the optimizer \code{\link[stats]{constrOptim}},
#'   including the input arguments (e.g., \code{k}, \code{options}, \code{A}, \code{V}, etc.).
#'   \itemize{
#'   \item If the Ab-representation was used, \code{par} provides the ML estimate for
#'   the probability vector \eqn{\theta}.
#'   \item If the V-representation was used, \code{par} provides the estimates for the
#'   (usually not identifiable) mixture weights \eqn{\alpha} that define the convex
#'   hull of the vertices in \eqn{V}, while \code{p} provides the ML estimates for
#'   the probability parameters. Because the weights must sum to one, the
#'   \eqn{\alpha}-parameter for the last row of the matrix \eqn{V} is dropped.
#'   If the unconstrained ML estimate is inside the convex hull, the mixture weights
#'   \eqn{\alpha} are not estimated and replaced by missings (\code{NA}).
#'   }
#'
#' @examples
#' # predicted linear order: p1 < p2 < p3 < .50
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
#' # probabilistic strategy:  A,A,A,B  [e1<e2<e3<e4<.50]
#' strat <- list(pattern = c(-1, -2, -3, 4),
#'               c = .5, ordered = TRUE, prior = c(1,1))
#' ml_binom(c(7,3,1, 19), 20, strategy = strat)[1:2]
#'
#'
#' # vertex representation (one prediction per row)
#' V <- matrix(c(
#'   # strict weak orders
#'   0, 1, 0, 1, 0, 1,  # a < b < c
#'   1, 0, 0, 1, 0, 1,  # b < a < c
#'   0, 1, 0, 1, 1, 0,  # a < c < b
#'   0, 1, 1, 0, 1, 0,  # c < a < b
#'   1, 0, 1, 0, 1, 0,  # c < b < a
#'   1, 0, 1, 0, 0, 1,  # b < c < a
#'
#'   0, 0, 0, 1, 0, 1,  # a ~ b < c
#'   0, 1, 0, 0, 1, 0,  # a ~ c < b
#'   1, 0, 1, 0, 0, 0,  # c ~ b < a
#'   0, 1, 0, 1, 0, 0,  # a < b ~ c
#'   1, 0, 0, 0, 0, 1,  # b < a ~ c
#'   0, 0, 1, 0, 1, 0,  # c < a ~ b
#'
#'   0, 0, 0, 0, 0, 0   # a ~ b ~ c
#'   ), byrow = TRUE, ncol = 6)
#' ml_multinom(k = c(4,1,5,  1,9,0,  7,2,1), n.fit = 1,
#'             options = c(3,3,3), V = V)
#' @export
ml_binom <- function(k, n, A, b, map, strategy, n.fit = 3, start,
                     progress = FALSE, ...){
  if (length(n) == 1) n <- rep(n, length(k))

  if (!missing(strategy) && !is.null(strategy)){
    pt <- strategy_to_Ab(strategy)
    A <- pt$A
    b <- pt$b
    map <- strategy$pattern
  } else {
    check_Abknprior(A, b, k, n)
    if (missing(map) || is.null(map))
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
ml_multinom <- function(k, options, A, b, V, n.fit = 3, start,
                        progress = FALSE, ...){

  t0 <- Sys.time()
  if (progress) cat("ML estimation starts at", format(t0), "using", n.fit, "starting values.\n",
                    "Current ML value: ")

  check_ko(k, options)
  n <- c(tapply(k, rep(1:length(options), options), sum))
  est_unconstr <- c(drop_fixed(k / rep(n, options), options))
  if (length(est_unconstr) != sum(options - 1))
    stop("ML estimation only implemented for a vetor k (not for data sets/matrices).")
  if (anyNA(est_unconstr)){
    inside_unconstr <- FALSE  # (n=0) required for prior sampling
  } else {
    ll_unconstr <- loglik_multinom(est_unconstr, k, options)
  }

  if(!missing(A) && !is.null(A)){
    ### Ab-representation
    check_Abokprior(A, b, options, k)
    tmp <- Ab_multinom(options, A, b, nonneg = TRUE)
    A <- tmp$A
    b <- tmp$b

    if (!anyNA(est_unconstr))
      inside_unconstr <- inside(est_unconstr, A = A, b = b)
    if (inside_unconstr){
      if (progress) cat("Unconstrained MLE within: A*x <= b;  -loglik:", ll_unconstr)
      oo <- list(par = est_unconstr,
                 value = ll_unconstr,
                 counts = NA,
                 convergence = 0,
                 message = "Unconstrained MLE satisfies constraint:  A*x <= b.",
                 outer.iterations = NA,
                 barrier.value = NA)
    } else {
      if (missing(start) || is.null(start) || !all(A %*% start < b))
        start <- find_inside(A, b, options = options)
      check_start(start, A, b, interior = FALSE)
      tryCatch (
        oo <- constrOptim(start, f = loglik_multinom, grad = grad_multinom,
                          k = k, options = options, ui = - A, ci = - b, ...),
        error = function(e) {
          print(e)
          cat("\n\n  (optimization failed with starting value:\n  ", start, ")\n")
        })
      cnt <- 1
      if (progress) cat("(",cnt,") ", round(oo$value, 1), "; ", sep = "")
      while (cnt < n.fit){
        start <- find_inside(A, b, random = TRUE)
        oo2 <- constrOptim(start, loglik_multinom, grad = grad_multinom,
                           k = k, options = options, ui = - A, ci = - b, ...)
        cnt <- cnt + 1
        if (progress) cat("(",cnt,") ",round(oo2$value, 1), "; ", sep = "")
        if (oo2$value < oo$value) oo <- oo2
      }
    }
    names(oo$par) <- colnames(A)
    oo$A <- A
    oo$b <- b


  } else {
    #### V-representation
    options <- check_V(V, options = options)
    S <- nrow(V)
    alpha_Ab <- Ab_multinom(options = S, nonneg = TRUE)

    if (!anyNA(est_unconstr))
      inside_unconstr <- inside_V(x = est_unconstr, V = V, return_glpk = TRUE)
    if (inside_unconstr$inside){
      if (progress) cat("Unconstrained MLE within convex hull of V;  -loglik:", ll_unconstr)
      oo <- list(par = rep(NA, S - 1),
                 value = ll_unconstr,
                 counts = NA,
                 convergence = 0,
                 message = "Unconstrained MLE is inside convex hull of V.",
                 outer.iterations = NA,
                 barrier.value = NA,
                 p = est_unconstr)
    } else {

      if(!missing(start) && !is.null(start) && (length(start) != (S - 1) || any(start < 0) || sum(start) >= 1))
        stop("'start' must be a vector of nonnegative mixture weights,\n",
             "   and defines the weights for the first ", S-1, " row vectors of 'V'.")
      else
        start <- rdirichlet(1, rep(3, S))[,-S]  # first starting value closer to center

      tryCatch (
        oo <- constrOptim(start, f = loglik_mixture, grad = grad_multinom_mixture,
                          k = k, options = options, V = V,
                          ui = - alpha_Ab$A, ci = - alpha_Ab$b, ...),
        error = function(e) {
          print(e)
          cat("\n\n  (optimization failed with starting value:\n  ", start, ")\n")
        })
      cnt <- 1
      if (progress) cat("(",cnt,") ", round(oo$value, 1), "; ", sep = "")
      while (cnt < n.fit){
        start <- rdirichlet(1, rep(1, S))[,-S]
        oo2 <- constrOptim(start, f = loglik_mixture, grad = grad_multinom_mixture,
                           k = k, options = options, V =V,
                           ui = - alpha_Ab$A, ci = - alpha_Ab$b, ...)
        cnt <- cnt + 1
        if (progress) cat("(",cnt,") ", round(oo2$value, 1), "; ", sep = "")
        if (oo2$value < oo$value) oo <- oo2
      }
    }
    names(oo$par) <- paste0("alpha[",1:(S-1), "]")
    if (!inside_unconstr$inside){
      alpha <- add_fixed(oo$par, options = S, sum = 1)
      oo$p <- c(t(V) %*% alpha)
    }
    names(oo$p) <- index_mult(options, fixed = FALSE)
    oo$V <- V
  }
  t1 <- Sys.time()
  if (progress) cat("\nFinished at", format(t1), " (difference:", format(t1-t0), ").\n")
  oo$k <- k
  oo$options <- options
  oo$n.fit <- n.fit
  oo$time <- t1-t0
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

# product binomial loglikelihood with luckiness [DEPRECATED]
#
# luckiness: -log(a(prob))   not normalized
# => complexity of order-constrained models can be compared
# default:    standard NML; luck=c(1,1)
# uniform BF: inverse of Jeffreys prior; luck=c(1.5, 1.5)
loglik_strategy <- function (error, k, n, strategy){
  pb <- error_to_prob(error, strategy)
  loglik_binom(pb, k, n)
}

loglik_binom <- function (p, k, n){
  lp1 <- log(p)
  lp1[k == 0] <- 0
  lp0 <- log(1 - p)
  lp0[k == n] <- 0
  ll <- sum(k * lp1, (n-k) * lp0)
  # try(ll <- sum(dbinom(x = k, size = n, prob = p, log = TRUE)), silent = TRUE)
  # if (is.null(ll) || !is.na(ll) && ll == - Inf && all(k <= n))
  #   ll <- MIN_LL
  if (is.na(ll) || ll == -Inf)
    warning("log-likelihood: NaN or -Inf: \n",
            "  Check whether model implies strictly positive probability for each k>0!")
  - ll
}

# log(k1)*p1 + log(k2)*p2 + ...
# log(k1)*(V[,1] %*% alpha') + ...   where: alpha' = c(alpha, 1-sum(alpha))
# NOT vectorized
loglik_mixture <- function (alpha, k, V, options){
  alpha <- c(alpha, 1 - sum(alpha))
  p <- c(t(V) %*% alpha)
  p_all <- add_fixed(c(p), options = options, sum = 1)
  lp <- log(p_all)
  lp[k == 0] <- 0
  ll <- sum(k * lp)
  if (is.na(ll) || ll == -Inf)
    warning("log-likelihood: NaN or -Inf: \n",
            "  Check whether model implies strictly positive probability for each k>0!")
  - ll
}

# loglik:    log(k1)*(V[,1] %*% c(alpha, 1-sum(alpha)) + ...
#            log(k1)* [sum V[s,1]*alpha[s] + V[S,1]*(1-sum(alpha))]

# p1 <- function(alpha){
#   S <- nrow(V)
#   alpha <- c(alpha, 1 - sum(alpha))
#   p <- t(V) %*% alpha   # = V[S,] + dp_dalpha %*% alpha[1:(S-1)]
#   1 - p[5] - p[6]
# }
# dp_dalpha[6,]
# -dpK_dalpha[3,]
# numDeriv::grad(p1, alpha[1:(S-1)])

# NOT vectorized
grad_multinom_mixture <- function(alpha, k, V, options){
  S <- nrow(V)
  alpha <- c(alpha, 1 - sum(alpha))
  p <- c(t(V) %*% alpha)   # = V[S,] + dp_dalpha %*% alpha[1:(S-1)]
  p_all <- add_fixed(c(p), options = options, sum = 1)
  idx_K <- cumsum(options)
  K <- k[idx_K]
  pK <- p_all[idx_K]

  # unlist: stability if V data.frame
  dp_dalpha <- t(V[1:(S-1),,drop=FALSE]) - unlist(V[S,])  # = d p_j(alpha) / d alpha (numderiv: ok)
  dpK_dalpha <- do.call("rbind", by(dp_dalpha, rep(1:length(options), options - 1),
                                    colSums))  # = d p_K(alpha) / d alpha (numderiv: ok)
  g <- c(k[-idx_K] / p) %*% dp_dalpha - (K /pK) %*% dpK_dalpha
  # check gradient:
  # nd <- numDeriv::grad(loglik_mixture, alpha[1:(S-1)], k=k, options=options, V=V)
  # print(rbind(g = - g, numderiv = nd, diff = g + nd))
  # ll <- sum(dmultinom(x = k, size = n, prob = p, log = TRUE))
  # if (!is.na(g) && ll == - Inf)
  # g[!is.na(g) & g == ] <- MIN_LL
  - g
}

loglik_multinom <- function (p, k, options, p_drop = FALSE){
  if (!p_drop){
    p_all <- add_fixed(p, options = options, sum = 1)
  } else {
    p_all <- p
  }
  lp <- log(p_all)

  if (is.null(dim(lp))){
    lp[k == 0] <- 0
    ll <- sum(k * lp)
  } else {
    lpt <- t(lp)
    lpt[k == 0] <- 0
    ll <- colSums(k * lpt)
  }
  # ll <- sum(k * lp)
  # tapply(p, oo, function(pp) dmultinom())
  # ll <- sum(dmultinom(x = k, size = n, prob = p, log = TRUE))
  # if (!is.na(ll) && ll == - Inf)
  #   ll <- MIN_LL
  if (anyNA(ll) || any(ll == -Inf))
    warning("log-likelihood: NaN or -Inf: \n",
            "  Check whether model implies strictly positive probability for each k>0!")
  - ll
}


grad_multinom <- function (p, k, options){
  p_all <- add_fixed(p, options = options, sum = 1)
  idx_K <- cumsum(options)
  K <- k[idx_K]
  pK <- p_all[idx_K]
  g <- k[-idx_K] / p - rep(K / pK, options - 1)
  # check: print(rbind(g = g, numderiv = numDeriv::grad(loglik_multinom, p, k=k, options=options)))
  # ll <- sum(dmultinom(x = k, size = n, prob = p, log = TRUE))
  # if (!is.na(g) && ll == - Inf)
  # g[!is.na(g) & g == ] <- MIN_LL
  - g
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

