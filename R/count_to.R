
#' Compute Bayes Factor Using Prior/Posterior Counts
#'
#' Computes the encompassing Bayes factor (and standard error) defined as the
#' ratio of posterior/prior samples that satisfy the order constraints (e.g., of a polytope).
#'
#' @param posterior a list with the entries \code{count} (number of samples within polytope)
#'     and \code{M} (total number of samples). Both can be vectors in the stepwise procedure.
#'     See \code{\link{count_binomial}}.
#' @param prior a list similar as \code{posterior} but based on prior samples.
#' @param beta prior parameters of beta distributions for estimating the standard error
#' @param samples number of samples from beta distributions used to approximate the standard error. The default is Jeffreys' prior.
#'
#' The encompassing model is the unconstrained baseline model that assumes a separate,
#' unconstrained probability for each observable frequency.
#'
#' The standard error of the (stepwise) encompassing Bayes factor is estimated by sampling
#' ratios from beta distributions, with parameters defined by the posterior/prior counts
#' (see Hoijtink, 2011; p. 211).
#'
#' @template ref_hoijtink2011
#' @template return_bf
#' @seealso \code{\link{count_binomial}}, \code{\link{count_multinomial}}
#' @examples
#' post  <- list(count = 1447, M = 5000)
#' prior <- list(count = 152, M = 5000)
#' count_to_bf(post, prior)
#' @export
count_to_bf <- function (posterior, prior, beta = c(.5, .5), samples = 3000){
  check_count(posterior)
  check_count(prior)
  const <- 0
  if (!is.null(posterior$const_map_0e))
    const <- posterior$const_map_0e
  est <- sum(log(posterior$count / posterior$M)) - sum(log(prior$count / prior$M)) + const

  bpost  <- sampling_beta(posterior$count, posterior$M, beta, samples)
  bprior <- sampling_beta(prior$count, prior$M, beta, samples)
  lbf_0e <- rowSums(log(bpost)) - rowSums(log(bprior)) + const
  lbf_e0 <- rowSums(log(bprior)) - rowSums(log(bpost)) + const

  bf <- matrix(
    c(exp(est),  exp(-est), est, - est,
      sd(exp(lbf_0e)), sd(exp(lbf_e0)), sd(lbf_0e), sd(lbf_e0)),
    4, 2,  dimnames = list(c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"),
                           c("BF", "SE")))
  # list("bf" = bf, "posterior" = posterior, "prior" = prior)
  bf
}


sampling_beta <- function (count, M, beta = c(.5, .5), samples = 10000){
  S <- length(count)
  if (length(M) == 1)
    M <- rep(M, S)
  probs <- matrix(NA, samples, S)
  for (s in 1:S){
    probs[,s]  = rbeta(samples, count[s] + beta[1], M[s] - count[s]  + beta[2])
  }
  probs
}
