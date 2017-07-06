#' Standard Error for Bayes Factor Estimate
#'
#' Computes the standard error of the (stepwise) encompassing Bayes factor by sampling
#' ratios from beta distributions, with parameters defined by the posterior/prior counts
#' (see Hoijtink, 2011; p. 211).
#'
#' @param npost number of posterior samples that satisfy order constraints
#'    (in the stepwise procedure: a vector of integers)
#' @param nprior number of prior samples that satisfy order constraints (see \code{npost})
#' @param Mpost total number of posterior samples that have been drawn to compute the BF.
#'    In the stepwise procedure: if the length of \code{Mpost} is one,
#'    it is assumed that the same number of draws was used for each step
#' @param Mprior total number of prior samples (cf. \code{Mpost})
#' @param prior prior parameters of the beta distribution of the posterior/prior counts.
#'    The default is Jeffreys' prior.
#' @param samples number of samples from beta distributions used to approximate the standard error
#' @template ref_hoijtink2011
#' @return a vector with standard errors for the Bayes factor (and log BF) of the order-constrained vs. encompassing (\code{bf_0e}) and the encompassing vs. order-constrained (\code{bf_e0}) model
#' @examples
#' npost <- 2453
#' nprior <- 153
#' M <- 10000
#' # Bayes factor in favor of order constraint (bf_0e):
#' npost / nprior
#'
#' # standard error of estimate:
#' se_bf(npost, nprior, M)
#' @export
se_bf <- function (npost, nprior, Mpost, Mprior = Mpost,
                   prior = c(.5, .5), samples = 5000){

  bpost  <- sampling_beta(npost,  Mpost,  prior, samples)
  bprior <- sampling_beta(nprior, Mprior, prior, samples)

  bf_0e <- apply(bpost,  1, prod) / apply(bprior, 1, prod)
  bf_e0 <- apply(bprior, 1, prod) / apply(bpost,  1, prod)

  c("bf_0e" = sd(bf_0e), "bf_e0" = sd(bf_e0),
    "log_bf_0e" = sd(log(bf_0e)), "log_bf_e0" = sd(log(bf_e0)))
}


sampling_beta <- function (n, M, prior = c(.5, .5), samples = 5000){
  S = length(n)
  if (length(M) == 1)
    M <- rep(M, S)
  probs <- matrix(NA, samples, S)
  for (s in 1:S){
    probs[,s]  = rbeta(samples, n[s] + prior[1], M[s] - n[s]  + prior[2])
  }
  probs
}
