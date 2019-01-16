#' Transform Bayes Factors to Posterior Model Probabilities
#'
#' Computes posterior model probabilities based on Bayes factors.
#'
#' @param ... one or more Bayes-factor objects for different models as returned
#'   by the functions \code{\link{bf_binom}}, \code{\link{bf_multinom}} and
#'   \code{\link{count_to_bf}} (i.e., a 3x4 matrix containing a row
#'   \code{"bf0u"} and a column \code{"bf"}). Note that the Bayes factors must
#'   have been computed for the same data and using the same prior (this is not
#'   checked internally).
#' @param prior a vector of prior model probabilities (default: uniform). The
#'   order must be identical to that of the Bayes factors supplied via
#'   \code{...}. If \code{include_unconstr=TRUE}, the unconstrained model is
#'   automatically added to the list of models (at the last position).
#' @param include_unconstr whether to include the unconstrained, encompassing
#'   model without inequality constraints (i.e., the saturated model).
#' @inheritParams count_to_bf
#'
#' @examples
#' # data: binomial frequencies in 4 conditions
#' n <- 100
#' k <- c(59, 54, 74)
#'
#' # Hypothesis 1: p1 < p2 < p3
#' A1 <- matrix(c(1, -1,  0,
#'                0,  1, -1), 2, 3, TRUE)
#' b1 <- c(0, 0)
#'
#' # Hypothesis 2: p1 < p2 and p1 < p3
#' A2 <- matrix(c(1, -1,  0,
#'                1,  0, -1), 2, 3, TRUE)
#' b2 <- c(0, 0)
#'
#' # get posterior probability for hypothesis
#' bf1 <- bf_binom(k, n, A = A1, b = b1)
#' bf2 <- bf_binom(k, n, A = A2, b = b2)
#' postprob(bf1, bf2,
#'          prior = c(bf1=1/3, bf2=1/3, unconstr=1/3))
#' @export
postprob <- function(..., prior, include_unconstr = TRUE){

  dots <- list(...)
  bfnames <- sapply(substitute(list(...))[-1], deparse)

  if(!all(sapply(dots, is.matrix)) ||
     !all(sapply(dots, nrow) == 3 &
          sapply(dots, ncol) == 4 &
     sapply(dots, function(x) "bf_0u" %in% rownames(x)) ) )
    stop("As input, Bayes-factor objects must be provided as those returned by bf_binom (i.e., 3x4 matrices).")

  bf0u <- sapply(dots, "[", "bf_0u", "bf")

  if (include_unconstr){
    bfnames <- c(bfnames, "unconstrained")
    bf0u <- c(bf0u, 1)
  }
  names(bf0u) <- bfnames

  if (missing(prior) || is.null(prior)){
    prior <- rep(1, length(bfnames))
  } else if (any(prior < 0) || length(prior) != length(bfnames)){
    stop("'prior' must be positive and have the same length as the number of BFs.")
  }
  prior <- prior / sum(prior)
  names(prior) <- bfnames

  postprob <- bf0u * prior / sum(bf0u * prior)
  cbind(prior = prior, posterior = postprob)
}

#### NOTE: different data structure required to get approximation error!
#
# => require counts (not bf!) for all models
# => use beta sampling
# => compute pp with SE
