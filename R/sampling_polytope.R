#' Posterior Sampling for Constrained Binomial Models
#'
#' Uses Gibbs sampling to draw posterior samples for
#' binomial models with linear inequality-constraints.
#'
#' @inheritParams count_binom
#' @param M number of posterior samples
#' @param burnin number of burnin samples that are discarded. can be chosen to be
#'     small if the maxmimum-likelihood estimate is used as a starting value.
#'
#' @details
#' Can be used to obtain posterior samples for binary random utility models that
#' assume a mixture over predefined preference orders/vertices that jointly define
#' a convex polytope via (A * x <= b).
#'
#' @return a matrix with posterior samples for the probabilities to choose Option B for each item type
#' @template ref_myung2005
#' @seealso \code{\link{count_binom}}
#' @examples
#' A <- matrix(c(1, 0, 0,   # x1 < .50
#'               1, 1, 1,   # x1+x2+x3 < 1
#'               0, 2, 2,   # 2*x2+2*x3 < 1
#'               0, -1, 0,  # x2 > .2
#'               0, 0, 1),  # x3 < .1
#'             ncol = 3, byrow = TRUE)
#' b <- c(.5, 1, 1, -.2, .1)
#' samp <- sampling_binom(c(5,12,7), c(20,20,20), A, b)
#' head(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0,.6))
#' @export
sampling_binom <- function (k, n, A, b, map = 1:ncol(A), prior = c(1, 1),
                            M = 5000, start, burnin = 10, progress = TRUE){
  if (length(n) == 1) n <- rep(n, length(k))
  if (length(prior) == 1) prior <- rep(prior, length(k))
  aggr <- map_k_to_A(k, n, A, map)
  k <- aggr$k
  n <- aggr$n
  check_Abknprior(A, b, k, n, prior)
  if (missing(start) || any(start < 0))
    start <-  ml_binom(k, n, A, b, map, n.fit = 1, maxit = 20)$par
  check_start(start, A, b, interior = TRUE)
  samples <- sampling_bin(k, n, A, b, prior, M, start, burnin, progress)
  colnames(samples) <- colnames(A)
  samples
}

#' Posterior Sampling for Constrained Multinomial Models
#'
#' Uses Gibbs sampling to draw posterior samples for
#' multinomial models with linear inequality-constraints.
#'
#' @inheritParams count_multinom
#' @inheritParams count_binom
#' @inheritParams sampling_binom
#' @examples
#' # binary and ternary choice:
#' #           (a1,a2   b1,b2,b3)
#' k       <- c(15,9,   5,2,17)
#' options <- c(2,      3)
#'
#' # columns:   (a1,  b1,b2)
#' A <- matrix(c(1, 0, 0,   # a1 < .20
#'               0, 2, 1,   # 2*b1+b2 < 1
#'               0, -1, 0,  # b1 > .2
#'               0, 0, 1),  # b2 < .4
#'             ncol = 3, byrow = TRUE)
#' b <- c(.2, 1, -.2, .4)
#' samp <- sampling_multinom(k, options, A, b)
#' head(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0, 1))
#' @export
sampling_multinom <- function (k, options, A, b, prior = rep(1, sum(options)),
                               M = 5000, start, burnin = 10, progress = TRUE){
  if (missing(start) || any(start < 0)){
    start <-  ml_multinom(k, options, A, b, n.fit = 1, maxit = 20)$par
  }
  check_start(start, A, b, interior = TRUE)
  if (missing(k) || (length(k) == 1 && k == 0))
    k <- rep(0, sum(options))
  if (length(prior) == 1) prior <- rep(prior, sum(options))
  check_Abokprior(A, b, options, k, prior)
  samples <- sampling_mult(k, options, A, b, prior, M, start, burnin, progress)
  colnames(samples) <- colnames(A)
  if (is.null(colnames(A)))
    colnames(samples) <- drop_fixed(index_mult(options), options)
  samples
}
