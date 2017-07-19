#' Posterior Samples for Ordered-Constrained Product-Binomial Models
#'
#' Uses Gibbs sampling to draw posterior samples for order-constrained models of binary choice.
#'
#' @inheritParams count_binomial
#' @param M number of posterior samples
#' @param start starting vectors. Must be in the convex polytope (that is, \code{A*start <= b})
#' @param burnin number of burnin samples that are discarded
#'
#' @details
#' Can be used to obtain posterior samples for binary random utility models that assume a mixture over predefined
#' preference orders/vertices that jointly define a convex polytope via (A * x <= b).
#'
#' @return a matrix with posterior samples for the probabilities to choose Option B for each item type
#' @template ref_myung2005
#' @seealso \code{\link{count_binomial}}
#' @examples
#' A <- matrix(c(1, 0, 0,   # x1 < .50
#'               1, 1, 1,   # x1+x2+x3 < 1
#'               0, 2, 2,   # 2*x2+2*x3 < 1
#'               0, -1, 0,  # x2 > .2
#'               0, 0, 1),  # x3 < .1
#'             ncol = 3, byrow = TRUE)
#' b <- c(.5, 1, 1, -.2, .1)
#' samp <- sampling_binomial(c(5,12,7), c(20,20,20), A, b)
#' head(samp)
#' colMeans(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0,.6))
#' @export
sampling_binomial <- function (k, n, A, b, prior = c(1, 1), M = 5000,
                               start, burnin = 10){
  check_Abknprior(A, b, k, n, prior)
  if (missing(start) || all(start < 0)){
    start <- -1  # TODO: find start value
  } else if (length(start) != ncol(A) || !all(A %*% start <= b)) {
    stop ("'start' must be in the convex polytope:  A*start <= b")
  }
  samples <- sampling_binomial_cpp(k, n, A, b, prior, M, start, burnin)
  colnames(samples) <- colnames(A)
  if (is.null(colnames(A)))
    colnames(samples) <- index_bin(k)
  samples
}

#' Posterior Samples for Ordered-Constrained Product-Multinomial Models
#'
#' Uses Gibbs sampling to draw posterior samples for order-constrained models of binary choice.
#'
#' @inheritParams count_multinomial
#' @inheritParams count_binomial
#' @inheritParams sampling_binomial
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
#' samp <- sampling_multinomial(k, options, A, b)
#' head(samp)
#' colMeans(samp)
#' apply(samp, 2, plot, type = "l", ylim = c(0, 1))
#' @export
sampling_multinomial <- function (k, options, A, b,
                                  prior = rep(1, sum(options)),
                                  M = 5000, start, burnin = 10){
  if (missing(start) || any(start < 0)){
    start <- -1  # TODO: find start value
  } else if (length(start) != ncol(A) || !all(A %*% start <= b)) {
    stop ("'start' must be in the convex polytope:  A*start <= b")
  }
  if (missing(k) || (length(k) == 1 && k == 0))
    k <- rep(0, sum(options))
  check_Abokprior(A, b, options, k, prior)
  tmp <- ineq_probabilities(options, A, b)
  A <- tmp$A
  b <- tmp$b
  samples <- sampling_multinomial_cpp(k, options, A, b, prior, M, start, burnin)
  colnames(samples) <- colnames(A)
  if (is.null(colnames(A)))
    colnames(samples) <- index_mult(options)
  samples
}
