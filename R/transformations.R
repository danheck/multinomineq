
#' Get NML Model Weights
#'
#' Computes the NML weights for the competing models
#' (comparable to posterior model probabilities).
#' @param nml vector or matrix (one row per person) with NML values
#' @export
wnml <- function(nml){
  UseMethod("wnml", nml)
}

#' @rdname wnml
#' @export
wnml.default <- function(nml){
  exp(-nml)/sum(exp(-nml))
}

#' @rdname wnml
#' @export
wnml.matrix <- function(nml){
  t(apply(nml, 1, wnml))
}

#' @rdname wnml
#' @export
wnml.data.frame <- function(nml){
  wnml(as.matrix(nml))
}


#' Get Posterior Model Probabilities
#'
#' Computes the posterior model probabilities based on the marginal likelihoods
#' and the prior model probabilities.
#' @param marginal vector of log-marginal probabilities
#' @param prior vector of prior model probabilities (default: uniform over models). Note that vector is normalized internally to sum to one.
#' @export
post_prob <- function(marginal, prior){
  UseMethod("post_prob", marginal)
}

#' @rdname post_prob
#' @export
post_prob.default <- function(marginal, prior){
  if (missing(prior))
    prior <- rep(1, length(marginal))
  exp(marginal)*prior/sum(exp(marginal) * prior)
}

#' @rdname post_prob
#' @export
post_prob.matrix <- function(marginal, prior){
  t(apply(marginal, 1, post_prob, prior = prior))
}

#' @rdname post_prob
#' @export
post_prob.data.frame <- function(marginal, prior){
  post_prob(as.matrix(marginal), prior = prior)
}


# strategy <- factor(modelNames[apply(NML, 1, which.min)], levels=modelNames)
