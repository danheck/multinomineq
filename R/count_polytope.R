#' Count Number of Prior/Posterior Samples in Polytope
#'
#' Fast C++ function to count the number of prior/posterior samples that fall into
#' the polytope defined via A*x <= b. Useful to compute the encompassing Bayes factor.
#'
#' @param k the number of Option B choices.
#'     The default \code{k=0} an \code{n=0} is equivalent to sampling from the prior.
#' @param n the number of choices per item type.
#' @inheritParams compute_bf
#' @return a list with the elements
#' \itemize{
#'     \item\code{integral}: estimated probability that samples are in polytope
#'     \item\code{count}: number of samples in polytope
#'     \item\code{M}: total number of samples
#' }
#' @examples
#' # linear order constraint:
#' # x1 < x2 < .... < x6 < .50
#' A <- matrix(c(1, -1, 0, 0, 0, 0,
#'               0, 1, -1, 0, 0, 0,
#'               0, 0, 1, -1, 0, 0,
#'               0, 0, 0, 1, -1, 0,
#'               0, 0, 0, 0, 1, -1,
#'               0, 0, 0, 0, 0, 1),
#'             ncol = 6, byrow = TRUE)
#' b <- c(0, 0, 0, 0, 0, .50)
#'
#' k <- c(0, 3, 2, 5, 3, 7)
#' n <- rep(10, 6)
#'
#' # check whether specific vector is in polytope:
#' A %*% c(.05, .1, .12, .16, .19, .23) <= b
#'
#' # count prior samples and compare to analytical result
#' prior <- count_polytope(A, b, M = 5e5)
#' prior
#' (.50)^6 / factorial(6)
#'
#' # count posterior samples and get Bayes factor
#' posterior <- count_polytope(A, b, k, n,
#'                       M=c(1e5, 2e4), steps = 2)
#' posterior$integral / prior$integral  # BF for constraints
#' count_to_bf(posterior, prior)
#' @export
count_polytope <- function (A, b, k = 0, n = 0, prior = c(1, 1),
                            M = 10000, steps, batch = 10000){
  if (length(k) == 1 && k == 0)
    k <- rep(0, ncol(A))
  if (length(n) == 1 && n == 0)
    n <- rep(0, ncol(A))
  check_knAbprior(k, n, A, b, prior)
  check_Mbatch(M, batch)

  if (missing(steps) || is.null(steps) || length(steps) == 0){
    cnt <- as.list(count_polytope_cpp(k, n, A, b, prior, M, batch))
  } else {
    check_stepsA(steps, A)
    cnt <- count_stepwise(k, n, A, b, prior, M, steps, batch)
  }
  cnt
}

#' Compute Bayes Factor Using Prior/Posterior Counts
#'
#' Computes the encompassing Bayes factor (and standard error) defined as the ratio of posterior/prior
#' samples that satisfy the order constraints (e.g., of a polytope).
#'
#' @param posterior a list with the entries \code{count} (number of samples within polytope)
#'     and \code{M} (total number of samples). Both can be vectors in the stepwise procedure.
#'     See \code{\link{count_polytope}}.
#' @param prior a list similar as \code{posterior} but based on prior samples.
#' @param beta prior parameters of beta distributions for estimating the standard error
#' @param samples number of samples from beta distributions used to approximate the standard error.
#'      The default is Jeffreys' prior.
#'
#' The standard error of the (stepwise) encompassing Bayes factor is estimated by sampling
#' ratios from beta distributions, with parameters defined by the posterior/prior counts
#' (see Hoijtink, 2011; p. 211).
#'
#' @template ref_hoijtink2011
#' @template return_bf
#' @examples
#' post  <- list(count = 1447, M = 5000)
#' prior <- list(count = 152, M = 5000)
#' count_to_bf(post, prior)
#' @export
count_to_bf <- function (posterior, prior, beta = c(.5, .5), samples = 3000){
  check_count(posterior)
  check_count(prior)
  est <- prod(posterior$count / posterior$M)  / prod(prior$count / prior$M)

  bpost  <- sampling_beta(posterior$count, posterior$M, beta, samples)
  bprior <- sampling_beta(prior$count, prior$M, beta, samples)

  bf_0e <- apply(bpost,  1, prod) / apply(bprior, 1, prod)
  bf_e0 <- apply(bprior, 1, prod) / apply(bpost,  1, prod)

  bf <- matrix(c(est,  1 / est, log(est), - log(est),
                 sd(bf_0e), sd(bf_e0), sd(log(bf_0e)), sd(log(bf_e0))),
               4, 2,  dimnames = list(c("bf_0e", "bf_e0", "log_bf_0e", "log_bf_e0"),
                                      c("BF", "SE")))
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

