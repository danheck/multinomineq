
#' Transform Bayes Factors to Posterior Model Probabilities
#'
#' Computes posterior model probabilities based on objects returned by
#' \code{\link{bf_binomial}}, \code{\link{bf_multinomial}} and \code{\link{count_to_bf}}.
#'
#' @param ... an arbitrary number of Bayes-factor objects with the elements \code{count} and \code{M} indicating how often posterior samples where within the admissible parameter space relative to the same encompassing model. Such lists are returned by \code{\link{count_binomial}} or \code{\link{count_multinomial}}.
#' @param prior prior model probabilities. Uniform by default.
#' @inheritParams count_to_bf
#' @export
postprob <- function(..., prior, beta = c(.5, .5), samples = 3000){
  dots <- list(...)
  if (missing(prior)){
    prior <- rep(1, length(dots))
  } else if (any(prior < 0) || length(prior) != length(dots)){
    stop("'prior' must be positive and have the same length as the number of BFs.")
  }
  prior <- prior / sum(prior)
  # TODO:
  # get counts for all models
  # use beta sampling
  # compute pp with SE
  NULL
}
