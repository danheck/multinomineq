
#' Compute Bayes Factor Using Prior/Posterior Counts
#'
#' Computes the encompassing Bayes factor (and standard error) defined as the
#' ratio of posterior/prior samples that satisfy the order constraints (e.g., of a polytope).
#'
#' @param posterior a vector (or matrix) with entries (or columns)
#'     \code{count} = number of posterior samples within polytope and \code{M} =
#'     total number of samples. See \code{\link{count_binom}}.
#' @param prior a vecotr or matrix similar as for \code{posterior}, but based on
#'     samples from the prior distribution.
#' @param beta prior parameters of beta distributions for estimating the standard error.
#' @param samples number of samples from beta distributions used to approximate
#'    the standard error. The default is Jeffreys' prior.
#'
#' The encompassing model is the unconstrained baseline model that assumes a separate,
#' unconstrained probability for each observable frequency. The Bayes factor is obtained
#' as the ratio of posterior/prior samples within an order-constrained subset of the
#' parameter space.
#'
#' The standard error of the (stepwise) encompassing Bayes factor is estimated by sampling
#' ratios from beta distributions, with parameters defined by the posterior/prior counts
#' (see Hoijtink, 2011; p. 211).
#'
#' @template ref_hoijtink2011
#' @template return_bf
#' @seealso \code{\link{count_binom}}, \code{\link{count_multinom}}
#' @examples
#' # vector input
#' post  <- c(count = 1447, M = 5000)
#' prior <- c(count = 152,  M = 10000)
#' count_to_bf(post, prior)
#'
#' # matrix input (due to nested stepwise procedure)
#' post <- cbind(count = c(139, 192), M = c(200, 1000))
#' count_to_bf(post, prior)
#' @export
count_to_bf <- function (posterior, prior, beta = c(.5, .5), samples = 3000){
  posterior <- check_count(posterior)
  prior <- check_count(prior)
  const <- 0
  if (!is.null(attr(posterior, "const_map_0e"))){
    const <- attr(posterior, "const_map_0e")
  }
  est <- sum(log(posterior[,1] / posterior[,2])) -
    sum(log(prior[,1] / prior[,2])) + const

  bpost  <- sampling_integral(posterior[,1], posterior[,2], log = TRUE, beta, samples)
  bprior <- sampling_integral(prior[,1], prior[,2], log = TRUE, beta, samples)
  lbf_0e <- bpost - bprior + const
  lbf_e0 <- bprior - bpost + const

  bf <- matrix(
    c(exp(est),  exp(-est), est, - est,
      sd(exp(lbf_0e)), sd(exp(lbf_e0)), sd(lbf_0e), sd(lbf_e0)),
    4, 2,  dimnames = list(c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"),
                           c("BF", "SE")))
  bf
}

precision_count <- function(count, M, log = TRUE,
                            beta = c(.5, .5), samples = 5000){
  s <- sampling_integral(count, M, log = TRUE, beta, samples)
  if (log){
    sd(s)
  } else {
    sd(exp(s))
  }
}

sampling_integral <- function(count, M, log = TRUE, beta = c(.5, .5), samples = 5000){
  S <- length(count)
  if (length(M) == 1)
    M <- rep(M, S)
  probs <- rep(0, samples)
  for (s in 1:S){
    probs <- probs + log(rbeta(samples, count[s] + beta[1],
                               M[s] - count[s]  + beta[2]))
  }
  if (log){
    probs
  } else {
    exp(probs)
  }
}

sampling_beta <- function (count, M, beta = c(.5, .5), samples = 5000){
  S <- length(count)
  if (length(M) == 1)
    M <- rep(M, S)
  probs <- matrix(NA, samples, S)
  for (s in 1:S){
    probs[,s]  = rbeta(samples, count[s] + beta[1], M[s] - count[s]  + beta[2])
  }
  probs
}
