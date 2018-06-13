
#' Get Posterior/NML Model Weights
#'
#' Computes the posterior model probabilities based on the log-marginal likelihoods/negative NML values.
#'
#' @param x vector or matrix of log-marginal probabilities or negative NML values (if matrix: one model per column)
#' @param prior vector of prior model probabilities (default: uniform over models). The vector is normalized internally to sum to one.
#' @examples
#' logmarginal <- c(-3.4, -2.0, -10.7)
#' model_weights(logmarginal)
#'
#' nml <- matrix(c(2.5, 3.1, 4.2,
#'                 1.4, 0.3, 8.2), nrow = 2, byrow = TRUE)
#' model_weights(-nml)
#' @export
model_weights <- function(x, prior){
  if (is.matrix(x) || is.data.frame(x)){
    w <- apply(x, 1, model_weights, prior = prior)
    return (t(w))
  } else {
    if (missing(prior) || is.null(prior))
      prior <- rep(1, length(x))
    w <- exp(x)*prior/sum(exp(x) * prior)
    return (w)
  }
}
