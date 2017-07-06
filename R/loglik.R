#' Choice Probabilities Implied by Strategy Predictions
#'
#' Returns a vector with probabilities of choosing Option B for a given set
#' of predictions and corresponding error probabilities.
#'
#' @param error parameter vector of error probabilities (depending on \code{prediction}).
#'   A single value for deterministic strategies and an
#'   order vector of probabilities for probabilistic strategies (cf. examples).
#' @inheritParams compute_cnml
#' @examples
#' # Predicted:  B,B,[guess],A
#' # (deterministic: with constant error = .10)
#' pred <- c(1, 1, 0, -1)
#' error_to_probB(error = .10, pred)
#'
#' # Predicted:  B,B,A,[guess]
#' # (probabilistic: with ordered errors e1<e2<e3)
#' pred2 <- c(1, 2, -3, 0)
#' error_to_probB(error = c(.10, .15,.20), pred2)
#' @export
error_to_probB <- function(error, prediction){
  # error per item type:
  if (!is.null(attr(prediction, "ordered"))){
    error_label <- rev(get_error_unique(prediction))
    error_per_type <- error[match(abs(prediction), error_label)]
  } else {
    error_per_type <- error
  }
  # reversed items and guessing:
  probB <- ifelse(prediction < 0, error_per_type,
                  ifelse(prediction > 0, 1 - error_per_type, .5))
  probB
}



# unique error indices/labels
get_error_unique <- function(prediction){
  pred <- prediction[prediction != 0]
  sort(unique(abs(pred)))
}

get_error_number <- function(prediction){
  length(get_error_unique(prediction))
}

# product binomial loglikelihood
loglik <- function (error, k, n, prediction){
  pb <- error_to_probB(error, prediction)
  ll <- sum(dbinom(x = k, size = n, prob = pb, log = TRUE))
  if (ll == - Inf && all(k <= n))
    ll <- MIN_LL
  ll
}

# count adherence/error rates
# used to get ML estimate : adherence/n
estimate_error <- function (k, n, prediction, luck = c(1,1), prob = TRUE){
  # reversed items:
  tmp <- k
  pred_B <- prediction > 0
  tmp[pred_B] <- n[pred_B] - k[pred_B]

  error_label <- get_error_unique(prediction)
  idx <- match(abs(prediction), error_label, nomatch = NA)
  cnt_predicted <- tapply(tmp, list(idx), sum) + luck[1] - 1
  if (prob){
    cnt_total <-  tapply(n, list(idx), sum) + sum(luck) - 2
    cnt_predicted <- cnt_predicted / cnt_total
  }
  unname(c(cnt_predicted), force = TRUE)
}

# move analytical estimate into interior of parameter space
adjust_error <- function(error, prediction, c = .50, bound = 1e-10){
  error_adj <- error
  o <- attr(prediction, "ordered")
  if (!is.null(o) && is.logical(o)){
    # simple linear order constraints
    if (o){
      error_adj <- adj_iterative(error, c, bound)
      error_adj <- sort(error_adj + runif(length(error), 0, bound/4))
    }
  } else if (!is.null(o)){
    #### complex order constraints (ui, ci)
    stop("complex order constraints (A, b) not yet implemented")
  } else if (length(error) > 0) {
    error_adj <- pmax(bound, pmin(c - bound, error))
  }
  error_adj
}

