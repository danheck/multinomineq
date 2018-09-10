
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
#' @param exact_prior optional: the exact prior probabability of the order constraints.
#'     For instance, \code{exact_prior=1/factorial(4)} if pi1<pi2<pi3<pi4 (and if  the prior is symmetric).
#'     If provided, \code{prior} is ignored.
#' @param log whether to return the log-Bayes factor instead of the Bayes factor
#' @param beta prior shape parameters of the beta distributions used for approximating the
#'     standard errors of the Bayes-factor estimates. The default is Jeffreys' prior.
#' @param samples number of samples from beta distributions used to compute
#'    standard errors.
#'
#'
#' The unconstrained (encompassing) model is the saturated baseline model that assumes a separate,
#' independent probability for each observable frequency. The Bayes factor is obtained
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
#'
#' # exact prior probability known
#' count_to_bf(posterior = c(count = 1447, M = 10000),
#'             exact_prior = 1/factorial(4))
#' @export
count_to_bf <- function (posterior, prior, exact_prior, log = FALSE,
                         beta = c(1/2, 1/2), samples = 3000){
  posterior <- check_count(posterior)
  const <- 0
  if (!is.null(attr(posterior, "const_map_0u"))){
    const <- attr(posterior, "const_map_0u")
  }
  s_post  <- sampling_proportion(posterior[,1], posterior[,2], log = TRUE,
                                 beta = beta, samples = samples)

  if (missing(exact_prior) || is.null(exact_prior)){
    prior <- check_count(prior)
    s_prior <- sampling_proportion(prior[,1], prior[,2], log = TRUE,
                                   beta = beta, samples = samples)

  } else{
    if (is.null(exact_prior) || !is.numeric(exact_prior) || length(exact_prior) != 1 ||
        exact_prior < 0 || exact_prior > 1)
      stop("'prior' must be a probability in the interval [0,1].")
    s_prior <- log(exact_prior)
    prior <- t(c(exact_prior, 1))
  }

  l_post <- sum(log(posterior[,1] / posterior[,2]))
  l_prior <- sum(log(prior[,1] / prior[,2]))
  est <- l_post - l_prior + const
  lbf_0u <- s_post - s_prior + const
  lbf_u0 <- s_prior - s_post + const

  # lbf_0n0 = log(f) - log(c) + log(1-c) - log(1-f)
  est_0n0 <- l_post - log1mexp(l_post) - l_prior + log1mexp(l_prior)
  lbf_0n0 <- s_post - s_prior + log1mexp(s_prior) - log1mexp(s_post)

  bf <- cbind("bf" = c("bf_0u" = est, "bf_u0" = -est, "bf_00'" = est_0n0),
              t(sapply(list(lbf_0u, lbf_u0, lbf_0n0),
                       summary_samples, exp = !log)))
  if(!log) bf[,"bf"] <- exp(bf[,"bf"])
  bf
}

# log(1 - exp(x))
# cf. https://cran.r-project.org/web/packages/Rmpfr/vignettes/log1mexp-note.pdf
log1mexp <- function(x, a0 = log(2)){
  ifelse(x <= a0,
         log(- expm1(x)),
         log1p(- exp(x)))
}

#' @importFrom stats quantile
summary_samples <- function(samples, probs = c(.05, .95), exp = FALSE){
  if(exp) samples <- exp(samples)
  qq <- quantile(samples, .05)
  c("se" = sd(samples),
    "ci" = quantile(samples, probs))
}

# check_bf <- function(bf){
#   is.matrix(bf) && dim(bf) == c(3,2) && identical(dimnames(bf), DIMNAMES_BF)
# }

precision_count <- function(count, M, log = TRUE,
                            beta = c(.5, .5), samples = 5000){
  s <- sampling_proportion(count, M, log = TRUE, beta, samples)
  if (log){
    sd(s)
  } else {
    sd(exp(s))
  }
}

sampling_proportion <- function(count, M, log = TRUE, beta = c(.5, .5), samples = 5000){
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

# sampling_beta <- function (count, M, beta = c(.5, .5), samples = 5000){
#   S <- length(count)
#   if (length(M) == 1)
#     M <- rep(M, S)
#   probs <- matrix(NA, samples, S)
#   for (s in 1:S){
#     probs[,s]  = rbeta(samples, count[s] + beta[1], M[s] - count[s]  + beta[2])
#   }
#   probs
# }
