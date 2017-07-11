#' Choice Probabilities Implied by Strategy Predictions
#'
#' Returns a vector with probabilities of choosing Option B for a predicted
#' pattern and corresponding error probabilities.
#'
#' @param error parameter vector of error probabilities (depending on \code{pattern}).
#'   A single value for deterministic strategies and an
#'   order vector of probabilities for probabilistic strategies (cf. examples).
#' @inheritParams compute_marginal
#' @examples
#' # Predicted:  B,B,[guess],A
#' # (deterministic: with constant error = .10)
#' s1 <- list(pattern = c(1, 1, 0, -1),
#'            ordered = FALSE, c = .1, prior = c(1,1))
#' error_to_prob(error = .10, s1)
#'
#' # Predicted:  B,B,A,[guess]
#' # (probabilistic: with ordered errors e1<e2<e3<.50)
#' s2 <- list(pattern = c(1, 2, -3, 0),
#'            ordered = TRUE, c = .5, prior = c(1,1))
#' error_to_prob(error = c(.10, .15,.20), s2)
#' @export
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
# luckiness: -log(a(theta))   not normalized
# => complexity of order-constrained models can be compared
# default:    standard NML; luck=c(1,1)
# uniform BF: inverse of Jeffreys prior; luck=c(1.5, 1.5)
loglik <- function (error, k, n, strategy){
  pb <- error_to_prob(error, strategy)
  ll <- sum(dbinom(x = k, size = n, prob = pb, log = TRUE))             # likelihood
  if (!is.null(error))
    ll <- ll + sum(dbeta(error, strategy$prior[1], strategy$prior[2], log = TRUE))  # luckiness
  if (!is.na(ll) && ll == - Inf && all(k <= n))
    ll <- MIN_LL
  ll
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

