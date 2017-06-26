#' Get Probabilities Predicted by Strategy
#'
#' Returns a vector with probabilities of choosing Option B for the different item types.
#'
#' @param par parameter vector (i.e., error probabilities)
#' @inheritParams compute_cnml
#' @examples
#' # Predicted: B, B, [guess], A
#' pred <- c(1, 1, 0, -1)
#' # error probability = .10
#' get_prob_B(par = .10, pred)
#'
#' # Predicted: B,B,A with errors: e1 < e2 < e3
#' pred2 <- c(1, 2, -3)
#' get_prob_B(par = c(.10, .15,.20), pred2)
#' @export
get_prob_B <- function(par, prediction){
  pb <- rep(NA, length(prediction))
  pb[prediction == 0] <- .50  # Guess = 0

  par_label <- get_par_unique(prediction)

  # replace pareter labels by current pareter values
  idxB <- match(prediction, par_label, nomatch = 0)
  pb[idxB != 0] <- 1 - par[idxB]
  idxA <- match(prediction, - par_label, nomatch = 0)
  pb[idxA != 0] <- par[idxA]
  pb
}

# unique parameter indices/labels
get_par_unique <- function(prediction){
  pred <- prediction[prediction != 0]
  sort(unique(abs(pred)))
}

get_par_number <- function(prediction){
  length(get_par_unique(prediction))
}

# product binomial loglikelihood
loglik <- function (par, k, n, prediction){
  pb <- get_prob_B(par, prediction)
  ll <- sum(dbinom(x = k, size = n, prob = pb, log = TRUE))
  if (ll == - Inf && all(k <= n))
    ll <- MIN_LL
  ll
}

# count adherence/error rates
# used to get ML estimate : adherence/n
estimate_par <- function (k, n, prediction, luck = c(1,1), prob = TRUE){
  # reversed items:
  tmp <- k
  pred_B <- prediction > 0
  tmp[pred_B] <- n[pred_B] - k[pred_B]

  par_label <- get_par_unique(prediction)
  idx <- match(abs(prediction), par_label, nomatch = NA)
  cnt_predicted <- tapply(tmp, list(idx), sum) + luck[1] - 1
  if (prob){
    cnt_total <-  tapply(n, list(idx), sum) + sum(luck) - 2
    cnt_predicted <- cnt_predicted / cnt_total
  }
  unname(c(cnt_predicted), force = TRUE)
}

# move analytical estimate into interior of parameter space
adjust_par <- function(par, prediction, c = .50, bound = 1e-10){
  par_adj <- par
  o <- attr(prediction, "ordered")
  if (!is.null(o) && is.logical(o)){
    # simple linear order constraints
    if (o){
      par_adj <- adj_iterative(par, c, bound)
      par_adj <- sort(par_adj + runif(length(par), 0, bound/4))
    }
  } else if (!is.null(o)){
    #### complex order constraints (ui, ci)
    stop("complex order constraints (ui, ci) not yet implemented")
  } else if (length(par) > 0) {
    par_adj <- pmax(bound, pmin(c - bound, par))
  }
  par_adj
}


