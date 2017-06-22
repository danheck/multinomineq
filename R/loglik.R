

# unique parameter labels
get_par_unique <- function(prediction){
  pred <- prediction[prediction != 0]
  sort(unique(abs(pred)))
}

get_par_number <- function(prediction){
  length(get_par_unique(prediction))
}

# model predictions: vector with probabilities to choose Option B
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


# product binomial loglikelihood
loglik <- function (par, b, n, prediction){
  pb <- get_prob_B(par, prediction)
  ll <- sum(dbinom(x = b, size = n, prob = pb, log = TRUE))
  #   ll <- sum(b*log(pb)+(n-b)*log(1-pb)) + lgamma(N+1) - lgamma(B+1) - lgamma(N-B+1) )
  if (ll == - Inf)
    ll <- MIN_LL
  ll
}


# get_adherence <- function(b, prediction){
#
#   par_label <- get_par_unique(prediction)
#   idx <- match(abs(prediction), par_label, nomatch = NA)
#   cnt_predicted <-  tapply(tmp, list(idx), sum)
#   unname(c(cnt_predicted), force = TRUE)
# }

estimate_par <- function (b, n, prediction, luck = c(1,1), prob = TRUE){
  # reversed items:
  tmp <- b
  pred_B <- prediction > 0
  tmp[pred_B] <- n[pred_B] - b[pred_B]

  par_label <- get_par_unique(prediction)
  idx <- match(abs(prediction), par_label, nomatch = NA)
  cnt_predicted <- tapply(tmp, list(idx), sum) + luck[1] - 1
  if (prob){
    cnt_total <-  tapply(n, list(idx), sum) + sum(luck) - 2
    cnt_predicted <- cnt_predicted / cnt_total
  }
  unname(c(cnt_predicted), force = TRUE)
}

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


