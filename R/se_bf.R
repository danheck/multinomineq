#' Standard Error for Bayes Factor Estimate
#'
#' Computes the standard error of the (stepwise) encompassing Bayes factor by sampling
#' ratios from beta distributions, with parameters defined by the posterior/prior counts
#' (see Hoijtink, 2011; p. 211).
#'
#' @param npost number of posterior samples that satisfy order constraints
#'    (in the stepwise procedure: a vector of integers)
#' @param npior number of prior samples that satisfy order constraints (see \code{npost})
#' @param M total number of posterior samples that have been drawn to compute the BF.
#'    In the stepwise procedure: if the length of \code{M} is one,
#'    it is assumed that the same number of draws was used for each step
#' @param prior prior parameters of the beta distribution of the posterior/prior counts.
#'    The default is Jeffreys' prior.
#' @param samples number of samples from beta distributions used to approximate the standard error
#' @template ref_hoijtink2011
#' @examples
#' npost <- 2453
#' nprior <- 153
#' M <- 10000
#' # Bayes factor in favor of order constraint:
#' npost / nprior
#'
#' # standard error of estimate:
#' se_bf(npost, nprior, M)
#' @export
se_bf <- function (npost, nprior, M,
                   prior = c(.5, .5), samples = 1000){
  S = length(npost)
  if (length(M) == 1) M <- rep(M, S)
  bpost = bprior = matrix(NA, samples, S)
  for (s in 1:S){
    bpost[,s] = rbeta(samples, npost[i] + prior[1], M[i] - npost[i] + prior[2])
    bprior[,s] = rbeta(samples, nprior[i] + prior[1], M[i] - nprior[i] + prior[2])
  }
  bfs <- apply(bpost, 1, prod) / apply(bprior, 1, prod)
  sd(bfs)
}
