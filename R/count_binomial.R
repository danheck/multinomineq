#' Number of Product-Binomial Prior/Posterior Samples in Polytope
#'
#' Draws prior/posterior samples for product-binomial data and counts how many samples are
#' inside the convex polytope defined by
#' (1) the inequalities A*x <= b or
#' (2) the convex hull over the vertices V.
#'
#' @inheritParams inside
#' @param k the number of Option B choices.
#'     The default \code{k=0} an \code{n=0} is equivalent to sampling from the prior.
#' @param n the number of choices per item type.
#' @param M number of posterior samples drawn from the encompassing model
#' @param batch size of the batches into which computations are split to reduce memory load
#' @param steps integer vector that indicates at which rows the matrix \code{A} is split for a stepwise computation of the Bayes factor (see details). In this case, \code{M} can be a vector with the number of samples drawn in each step from the (partially) order-constrained models using Gibbs sampling
#' @param prior a vector with two positive numbers defining the shape parameters of the beta prior distributions for each binomial rate parameter.
#' @param start only if \code{steps} is defined: a vector with starting values
#'     within the polytope. If \code{steps = -1}, a random starting value is used
#'     (but a poitn within the polytope might not be found if the order constraints are strong).
#'
#' @details
#' Useful to compute the encompassing Bayes factor for testing order constraints (see \code{\link{bf_binomial}}; Hojtink, 2011).
#'
#' The stepwise computation of the Bayes factor proceeds as follows:
#' If the steps are defined as \code{steps=c(5,10)}, the BF is computed in three steps by comparing:
#' (1) Unconstrained model vs. inequalities in \code{A[1:5,]};
#' (2) use posterior based on inequalities in \code{A[1:5,]} and check inequalities \code{A[6:10,]};
#' (3) sample from A[1:10,] and check inequalities in \code{A[11:nrow(A),]} (i.e., all inequalities).
#'
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
#' prior <- count_binomial(0, 0, A, b, M = 5e5)
#' prior
#' (.50)^6 / factorial(6)
#'
#' # count posterior samples and get Bayes factor
#' posterior <- count_binomial(k, n, A, b,
#'                             M=c(1e5, 2e4), steps = 2)
#' posterior$integral / prior$integral  # Bayes factor
#' count_to_bf(posterior, prior)
#' @template ref_hoijtink2011
#' @template ref_fukuda2004
#' @importFrom Rglpk Rglpk_solve_LP
#' @export
count_binomial <- function (k = 0, n = 0, A, b, V, prior = c(1, 1),
                            M = 10000, steps, batch = 10000, start = -1){

  if (missing(V) || is.null(V)){
    if (length(k) == 1 && k == 0)
      k <- rep(0, ncol(A))
    if (length(n) == 1 && n == 0)
      n <- rep(0, ncol(A))
    check_Abknprior(A, b, k, n, prior)
    check_Mbatch(M, batch)

    if (missing(steps) || is.null(steps) || length(steps) == 0){
      cnt <- as.list(count_binomial_cpp(k, n, A, b, prior, M, batch))
    } else {
      check_stepsA(steps, A)
      cnt <- count_stepwise(k, n, A, b, prior, M, steps, batch, start)
    }

  } else {
    if (length(k) == 1 && k == 0)
      k <- rep(0, ncol(V))
    if (length(n) == 1 && n == 0)
      n <- rep(0, ncol(V))
    count <- 0
    m <- M
    a <- c(rbind(k + prior[1], n - k + prior[2]))
    while (m > 0 ){
      X <- rpdirichlet_free(m, a, rep(2, ncol(V)))
      count <- count + sum(inside_V(X, V))
      m <- m - batch
    }
    cnt <- list("integral" = count/M, "count" = count, "M" = M)
  }
  cnt
}

# count_V <- function(X, V){
#   apply(X, 1, inside_V, V = V)
# }


#' Compute Bayes Factor Using Prior/Posterior Counts
#'
#' Computes the encompassing Bayes factor (and standard error) defined as the ratio of posterior/prior
#' samples that satisfy the order constraints (e.g., of a polytope).
#'
#' @param posterior a list with the entries \code{count} (number of samples within polytope)
#'     and \code{M} (total number of samples). Both can be vectors in the stepwise procedure.
#'     See \code{\link{count_binomial}}.
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
#' @seealso \code{\link{count_binomial}}, \code{\link{count_multinomial}}
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

